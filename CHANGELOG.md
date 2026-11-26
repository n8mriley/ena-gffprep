# Changelog

## Unreleased — User documentation

- Add release-pinned installation, offline quick start and guide navigation.
- Add compatibility reporting guidance, contributor instructions and a bug report template.

## Unreleased — Submission metadata and documentation

- Require BioProject for preparation via CLI, configuration, batch or direct API; always write STUDY in generated manifests. Audit remains exempt.
- Record BioProject in preparation metadata; validate syntax locally, not study ownership.
- Clarify primary ab initio/deep-learning scope and experimental evidence-based compatibility; streamline guides and document module responsibilities.

## Unreleased — ENA-only metadata

- Replace NCBI assembly metadata retrieval and parsing with ENA sequence-report TSV.
- Use versioned GCA accessions and separate ENA caches; remove RefSeq output selection.
- Rename `--no-ncbi-lengths` to `--no-ena-lengths` and report keys to ENA equivalents.
- Preserve local length inputs and manual organelle exclusions.

## Unreleased — Manual organelle exclusions

- Accept an optional comma-separated sequence list after `--remove-organelles`;
  retain automatic identification when the flag has no value (now ENA).
- Match exact input IDs, remove features and sequence declarations before
  normalization/tag allocation, and retain relationship and tag safeguards.
- Record manual provenance and reject unknown, duplicate or empty IDs.
- Support manual lists in JSON configuration and batch TSV; clarify that missing
  sequence lengths can still require remote metadata unless lookup is disabled.
- Add eight regression tests (210 total).

## 0.8.0

- Replace source-selected transformations with one source-independent graph normalizer; keep parsing and provenance separate.
- Add input-format, normalization and CDS-grouping policies; accept earlier source-format hints for compatibility.
- Normalize unknown annotators using the same rules, preserve explicit distinct protein products, and reject ambiguous/conflicting merges.
- Expand shared leaves before CDS grouping; retain identifier provenance and report applied structural rules.
- Preserve valid prior annotation outputs while accepting compatible existing CDS IDs without annotator-specific naming requirements.
- Keep organelle filtering before structural rewrites and tag allocation.

## 0.7.0

- Reorganize documentation around the general workflow, consolidate repeated guidance and add multi-annotator examples.

- Fall back to exact-assembly NCBI sequence lengths when local lengths are missing; share metadata/cache with organelle filtering.
- Add offline sequence metadata, shared cache/refresh options and a switch to disable length fallback; record lookup hashes and alias mappings.
- Detect supported source/dialect signatures from file contents, separately from pipeline provenance.
- Add Tiberius, mixed Tiberius/TransDecoder, Carbon-A, OrionGeno GTF, Helixer and BRAKER-compatible normalization.
- Preserve isoforms, source fields, coordinates, scores, CDS phases and existing evidence; reject ambiguous CDS product grouping.
- Add optional per-source inference JSON mapping and explicit prediction-method policy; direct ANNEVO/Helixer/Tiberius/Carbon-A/OrionGeno sources use ab initio provenance while other components retain existing evidence unless explicitly configured.
- Normalize non-CDS phases and child transcript identifiers with original values retained.
- Preserve ANNEVO conversion behavior; write deterministic gzip headers.
- Document format limits and separate local validation from Webin acceptance.

## 0.6.1 — Documentation review

- Focus validation documentation and verification metadata on the current test suite and supported checks.
- Remove nested verification history and dataset-specific test narratives.
- Clarify custom descriptions, unversioned ANNEVO headers, organelle exclusions and locus-tag reservations.
- Update the software version; annotation preparation behaviour is unchanged.

## 0.6.0 — Optional organelle exclusion

- Add --remove-organelles with exact-assembly NCBI lookup, validated caching, refresh and offline metadata input.
- Exclude organelle features before conversion and tag allocation; preserve reserved tags and check retained references.
- Add exclusion reports and CLI/config/batch support; keep audit mode descriptive.
- Format default descriptions with Annotator, Version and optional Model; preserve custom descriptions.

## 0.5.6 — Predictor and model metadata

- Select ANNEVO processing from the predictor in automatic mode and recognise unversioned headers.
- Add --predictor-model and accession.predictor.version[.model] default naming.
- Accept supported spacing and case variants in ANNEVO length comments.

## 0.5.5 — Documentation and report naming

- Clarify supported scope and rename the feature summary to feature_types_in_supported_subset.

## 0.5.4 — GFF-focused distribution

- Focus scripts and examples on GFF preparation and diagnostic review.

## 0.5.3 — Locus tags

- Preserve alphanumeric suffixes, reserve historical assignments and exclude repeat_region from propagation.

## 0.5.2 — Code documentation

- Add module, function and class documentation.

## 0.5.1 — Optional taxonomy

- Make taxon metadata optional while checking existing and explicit values for conflicts.

## 0.5.0 — Streamlined workflow

- Add package naming, manifests, reports, predictor detection, configuration and batch preparation.
