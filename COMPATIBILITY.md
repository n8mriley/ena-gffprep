# Compatibility and validation evidence

The primary scope is preparation of deep-learning and ab initio nuclear genome
annotations. Recognition of a source label is not a guarantee of compatibility
with every version, model or output mode.

## Evidence levels

| Level | Meaning |
| --- | --- |
| Synthetic regression | Automated checks of specific parsing, normalization or workflow behavior. |
| Real-data preparation | A real annotation was prepared successfully with recorded settings. |
| Feature concordance | Retained biological features were compared with the input. |
| Webin validation | The exact prepared package passed the intended Webin validation workflow. |
| Archive acceptance | ENA accepted the submission and assigned an accession. |

These are distinct outcomes, not interchangeable claims. BRAKER and other
evidence-based or mixed annotations remain experimental; synthetic coverage does
not establish end-to-end compatibility. See [scope](README.md#scope-and-maturity).

## Record a real-data test

Use one record per assembly, annotation and preparation configuration. Record
unknown or untested values explicitly; do not infer success from missing logs.
Only publish data and identifiers cleared for public release.

| Field | Value to record |
| --- | --- |
| Species and assembly | Scientific name and exact versioned GCA accession |
| Predictor | Name, version, model and ab initio/evidence-assisted mode |
| Input | Format, checksum and public source if available |
| Preparation | Toolkit tag/commit, command or config, metadata report checksum |
| Filtering | Excluded sequence IDs and gene/feature counts |
| Gene counts | Input and retained counts, accounting for exclusions |
| Feature concordance | Compared feature types, coordinates, strands, CDS phases and product relationships; counts of differences |
| Webin | Exact build, validation date, pass/fail/not run, warning summary |
| Archive | Annotation accession if accepted; otherwise not submitted or pending |

Compare retained features after accounting for documented ID remapping, shared
leaf expansion and normalization. A raw line count or byte comparison is not a
biological equivalence test. Check transcript membership and distinct products,
not only coordinate sets. Investigate unexpected differences before submission.
Retain preparation reports and validation logs with the test record, removing
credentials and confidential information before publication.

## Published real-data results

No complete real-data validation records are included in this guide yet. This
is a reporting template, not a benchmark or a claim of archive acceptance.
