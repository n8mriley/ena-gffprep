# Release verification — 0.6.0

## Current checks

All 120 automated test methods pass using Python standard-library unittest.
Tests never contact NCBI or Webin. Run:

```bash
python3 -m unittest discover -s tests -v
```

Existing GFF conversion, inference, metadata, locus-tag and output-publication
checks remain. New tests cover mitochondrial/plastid removal, exact versioned
alias matching, no-organelle/no-annotation cases, unmatched and ambiguous IDs,
unknown/missing classifications, incomplete JSON/JSONL, length conflicts,
cross-boundary references, all-features removal, reserved tags, no-flag parity,
cache reuse/refresh, pagination, transient/permanent network failures, corrupt
cache, config/batch handling, package failure and exclusion sidecars.

The latest exact Annotator / Version / Model description and unchanged CDS
inference are also covered.

## Real-data integration

The supplied Abraxas sylvata ANNEVO 2.2.1 annotation was prepared using the live
NCBI sequence report for GCA_964187985.1, model Insecta and a local test prefix.
The metadata contained 60 sequence records and classified OZ076621.1
(scaffold_MT_1) as Mitochondrion. All annotated sequences matched.

- Input: 263,764 feature rows and 16,528 genes.
- Removed: 50 rows (7 genes, 7 mRNAs, 18 exons, 18 CDS intervals).
- Retained: 263,714 feature rows and 16,521 genes.
- Output sequence-region declarations omit OZ076621.1.
- All retained feature rows exactly match the earlier unfiltered output, including tags.
- A second offline preparation with the previous cumulative map reserves all seven excluded gene tags.

The resulting file passed local gff3tools validation using the supplied
webin-cli-10.0.0-beta.jar: exit 0, no errors, 6,945 EXON_LENGTH warnings and
32 CDS_INTRON_LENGTH warnings (6,977 total). No validation rules were disabled.
This is bundled local validation, not authenticated Webin validation, submission
or archive-processing acceptance.

Earlier unfiltered user-provided validation results and older integration
results remain historical evidence, not new results for every release.

## Limits

The filter uses metadata, not nucleotide comparison or organelle discovery.
Unlabelled organelles and misclassified sequences can remain. Primary-assembly
unplaced/unlocalized scaffolds lacking a compartment are explicitly reported as
unclassified and retained with a warning. Other uncertain classifications fail.
Offline JSONL completeness depends on the supplied file; all annotated sequence
IDs must match. No registered prefix/account/assembly-availability checks or
network submission are performed by preparation. See ORGANELLES.md for policy.
