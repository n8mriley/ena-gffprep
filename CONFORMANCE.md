# Submission conformance (development)

Target: INSDC GFF3 v0.51 working draft, Submission [S] profile, based on SO GFF3 v1.26. This is a preparation tool, not a complete conformance validator. Exchange and Delivery requirements are not imposed on submissions.

## Local checks

In addition to the existing syntax, graph, coordinate, locus-tag and manifest checks, the development update checks taxonomy URI providers, Dbxref syntax and registered database tags, explicit amino acid names, multipart partial flags in biological order, and CDS-child exceptional stop codons. The amino acid list covers the 20 standard residues plus selenocysteine and pyrrolysine; other representations require review. Unknown pragmas/directives are omitted from output and recorded in the report; ordinary comments follow existing policies.

A single `#!transl_table` pragma accepts comma-separated compartment/code pairs, following the working draft directive syntax. Existing CDS-level values take precedence. For example, `#!transl_table primary:1,mitochondrion:2` is expanded into per-CDS attributes using exact ENA sequence aliases and lengths. These numbers are illustrative, not universal genetic codes. The draft does not enumerate compartment keys; this implementation uses `primary` for explicitly nuclear primary-assembly records and lower-case ENA organelle types (colon replaced by underscore). Unclassified scaffolds, contradictory metadata and missing compartment codes require explicit CDS values. A primary-only pragma without available metadata retains the previous primary-annotation assumption. Explicit CDS values need no compartment lookup. Without a pragma, the existing configured/default fallback remains unchanged; this is not taxonomic verification.

Use `--sequence-metadata` for offline classification, or enabled ENA length lookup for remote/cache resolution. Metadata already loaded for lengths or filtering is reused. Reports include input declarations, resolved compartments, metadata provenance, code counts per sequence and counts by assignment source. No sequence translation or BioSample genetic-code verification is performed.

## Offline vocabulary snapshot

Feature mappings are generated locally from a user-acquired gff3tools 5.5.1 JAR; Dbxref tags are supplied locally from the INSDC registry. No registry or extracted mapping is bundled. The setup tool checks the pinned JAR and OBO checksums and records the tag-list checksum. Tags are case-sensitive; identifiers are not remotely verified. The local tag list is user-maintained and its completeness is not verified. Pseudogene subtype normalization and explicit structural checks remain unchanged.

Selected conventional gene relationships are checked, including CDS parents and mRNA/gene relationships. This is not exhaustive ontology part_of validation. Existing prediction-method settings remain available; no new inference semantics are introduced by this update.

## Specialist feature semantics

Table 5 of the supplied working draft permits ribosomal slippage exceptions to CDS splice coordinates and recommends ordinary feature types with pseudogene attributes. The preparation engine now preserves explicitly declared one-base CDS interval overlaps for `exception=ribosomal_slippage`, on either strand. Duplicate intervals and larger overlaps remain refused. No slippage is inferred and no reading frame is repaired; exception groups remain excluded from ordinary phase-continuity review.

The supported `pseudogenic_CDS`, `pseudogenic_tRNA`, `pseudogenic_rRNA` and `pseudogenic_transcript` names are converted to their ordinary feature types with `pseudo=true` before phase checking and structural normalization. Coordinates and CDS phases are retained. Conflicting pseudo flags are refused. Repeated recoded-codon and anticodon IDs must agree on their amino acid attributes. Split codons remain subject to the three-base and parent-interval checks; these checks do not establish translation-frame alignment.

Additional RNA and exon types require an explicit gene owner under the local genome-annotation policy. This is a conservative preparation policy, not an assertion that the Submission profile universally mandates a gene parent. The SO graph is generated locally.

The origin-spanning section requires authoritative sequence-record topology. ENA sequence reports used here do not establish topology; neither an organelle label nor `Is_circular=true` alone is sufficient. Virtual coordinates therefore remain refused. Explicit `part`, remote locations, trans-splicing, circular RNA and embedded peptide FASTA also remain refused pending safe handling throughout parsing, normalization, filtering and output.

## Remaining gaps

- Complete parent relationships and broader gene ownership remain outstanding. The local feature vocabulary is a fixed subset; a valid term absent from it is refused for specialist review, not declared biologically invalid.
- Partial checks cover repeated logical IDs, not every possible exon/parent partial relationship. Sequence-based biological partiality and codon alignment require further checks.
- Circular virtual coordinates, trans-splicing, remote locations, explicit interval order and embedded peptide FASTA remain unsupported. No biological repair is attempted.
- Existing automatic predictor inference policies remain unchanged; explicitly choose prediction method and evidence policy for assisted annotations.
- Full authenticated Webin validation and real-data compatibility remain external requirements.

## External validation

Use a compatible Webin build. The standalone helper supports optional reference sequence, JSON metrics, JVM heap and processor settings, and a JSON summary with validator identity, JAR checksum and dynamically counted rule identifiers. Standalone defaults are not identical to Webin rules; no Webin-equivalent profile is claimed. See VALIDATION.md.

Sources: [SO GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md), [ENA feature mapping](https://github.com/enasequence/gff3tools/blob/main/src/main/resources/feature-mapping.tsv), [INSDC Dbxref vocabulary](https://www.insdc.org/submitting-standards/dbxref-qualifier-vocabulary/). The supplied working draft is not redistributed.

### Stage 5: pinned Sequence Ontology graph

`so_snapshot.json` is generated locally from the user-acquired SO OBO release 2026-08-07. It is not distributed. See README and THIRD_PARTY_NOTICES.md for setup and provenance. No runtime download or extra Python dependency is required after setup.

The working GFF3 specification column 3 requires a sequence feature or a term
related to SO:0000110 through is_a or part_of, and separately references the INSDC
feature mapping. We check both layers: obsolete and non-feature terms are refused;
valid but unmapped terms still require specialist review. Existing validator
mapping aliases remain supported (including legacy rRNA names); these are not
claims of canonical SO naming. The report records ontology version, checksum,
gene/transcript subtype names and parent relationship evidence counts per logical
feature. Positive evidence uses asserted part_of edges inherited through is_a;
absence is not proof of invalidity. Existing explicit structural checks remain
in force. This is not an OWL reasoner or an exhaustive INSDC mapping. New gene
subtypes are not automatically admitted or rewritten. No coordinate or phase
changes are introduced. The INSDC Minimal Specifications archive (commit
2a8668680625025178a8a38020ee426a4cc6070c) remains the reference for subsequent
biological and topology checks, not a substitute for the feature mapping.
