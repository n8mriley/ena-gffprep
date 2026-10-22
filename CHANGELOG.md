# Changelog

## 0.7.0

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
- Remove the packaged installation and repository administration guide.
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
