# Locus tags

A locus tag identifies a gene. Its supported transcript/RNA, CDS, exon, intron
and UTR components share the same tag, including alternative transcripts.
Allocation spans all sequences and gene types in an annotation. `repeat_region`
is excluded; an input repeat_region carrying a tag must be corrected explicitly.

## Prefix registration and tag syntax

Register the prefix through production ENA Webin's study registration workflow,
selecting “Will you provide functional genome annotation?”. Allow approximately
24 hours for a newly registered prefix to become usable. The converter checks
syntax, not registration, availability or ownership.

| Element | Converter rule |
| --- | --- |
| Prefix | `[A-Z][A-Z0-9]{2,11}`: 3–12 uppercase alphanumeric characters, starting with a letter. |
| Separator | One underscore between prefix and suffix. |
| Generated suffix | Positive integer, zero-padded to a minimum of six digits by default. |
| Existing suffix | ASCII alphanumeric, preserved exactly; examples: `WSITEST_A0001`, `ABC_r1112`. |
| Rejected forms | Missing separator, punctuation/decimal versions such as `ABC_00001.1`, nonpositive numeric suffixes and ambiguous differently padded versions of the same number. |

Use the registered namespace consistently across the assembly. Different
chromosomes or RNA types do not need separate prefixes. Five letters or five
suffix digits are example conventions, not universal requirements. A durable
prefix avoids dependence on a species or strain name that might change.

Use `--width 5` if five-digit generated suffixes are preferred. Width is a
minimum; numbers expand instead of truncating or wrapping. The README lists
[allocation options](README.md#locus-tags), including start and step.

## First allocation and updates

Prepare the whole assembly annotation together. First-time numbering follows
deterministic sequence-ID and coordinate ordering, not necessarily input row
order. Tags supplied in input attributes are checked for compatibility.

Keep the generated `locus_tags.tsv` permanently. On each update, add:

```bash
--locus-tag-map previous-package/locus_tags.tsv
```

This is an option on the normal preparation command, not a standalone command.
The map has `gene_id` and `locus_tag` columns; generated maps also contain
`status` (`active` or `reserved`). Use the latest complete cumulative map:

- Existing gene IDs retain their exact assignments, regardless of coordinate order.
- Missing genes' tags remain reserved in the new map and are not automatically reused.
- New numeric tags are allocated above the active/reserved numeric high-water mark,
  subject to start/step. Alphanumeric suffixes are reserved without raising that maximum.
- Input tags and the map must agree; one tag cannot belong to different gene IDs.

Matching uses source gene IDs, not coordinates or sequence similarity. Renamed,
split or merged genes require curated reconciliation. Omitting the historical
map removes the tool's ability to protect assignments absent from the input.

DDBJ permits insertion into unused numeric gaps. This tool does not fill gaps
automatically because they may represent retired tags. An independently verified
unused assignment can be supplied explicitly in the map while retaining all
historical reservations. Independent batch jobs with the same prefix are refused
to prevent colliding allocations.

## Filtering and feature scope

Organelle filtering happens before new tags are assigned. Excluded genes receive
no new tags; historical map entries and resolvable tags on excluded features
remain reserved. Keep the cumulative map across filtered and unfiltered runs.
See [ORGANELLES.md](ORGANELLES.md).

RNA and specialist components still need a supported gene hierarchy. Recognition
for tag propagation is not a guarantee that the intended Webin build accepts a
feature term. The preparation report records the allocation policy; external
validation remains necessary.

References:

- [ENA study and locus-tag registration](https://ena-docs.readthedocs.io/en/latest/submit/study/interactive.html)
- [DDBJ locus-tag conventions](https://www.ddbj.nig.ac.jp/ddbj/locus_tag-e.html)
