# Locus tags

A locus tag identifies a gene, not an individual exon, transcript or CDS interval.
All supported components of that gene inherit the same tag. Protein-coding and
non-coding RNA genes share one allocation sequence across every chromosome,
plasmid and other sequence in the input annotation. Alternative transcripts do
not receive separate locus tags.

## Registration and syntax

Register your prefix through the production ENA Webin study portal, selecting
“Will you provide functional genome annotation?”. ENA documents a delay of
approximately 24 hours before a new prefix becomes usable. Registration is
external: the scripts cannot verify ownership, availability or registration.
DDBJ's BioSample registration instructions describe its own submission route;
use the ENA registration route for this workflow.

Prefixes must match `[A-Z][A-Z0-9]{2,11}`: 3–12 uppercase alphanumeric characters,
starting with a letter. Digits are allowed after the initial letter, e.g. A1C.
Use the same registered prefix for all sequences in one assembly. Do not create
separate prefixes for chromosomes or RNA types. Choose a durable prefix rather
than encoding a species or strain name that may change.

The underscore between prefix and suffix is required. The tool accepts existing
suffixes matching `[A-Za-z0-9]+`, including `WSITEST_A0001` and the DDBJ example
`ABC_r1112`. Numeric generation remains the default:

```bash
python ena_gffprep.py \
  --input annotation.gff3 \
  --assembly-accession GCA_xxxxxxxxx.x \
  --locus-tag-prefix WSITEST \
  --predictor YOUR_PREDICTOR \
  --predictor-version YOUR_VERSION \
  --width 5
```

Use the actual predictor and version. A recognised versioned ANNEVO header
can supply them automatically; an unversioned header does not provide a version. This example generates `WSITEST_00001`, `WSITEST_00002`, etc.
The default numeric width is six digits. Width is a minimum: allocation
expands if the number outgrows it rather than truncating or wrapping.

Do not treat five letters or five suffix digits as universal requirements: those
are example conventions, not the documented prefix rules. Decimal versions such
as `WSITEST_00001.1`, missing separators and punctuated suffixes are rejected.
Existing numeric suffixes must be positive; differently padded representations
of the same number remain rejected as ambiguous by this tool.

## Updating an annotation

Keep the complete `locus_tags.tsv` permanently, including reserved entries, and
supply it on every update:

```bash
python ena_gffprep.py \
  --input updated.gff3 \
  --assembly-accession GCA_xxxxxxxxx.x \
  --locus-tag-prefix WSITEST \
  --predictor YOUR_PREDICTOR \
  --predictor-version YOUR_VERSION \
  --locus-tag-map previous-package/locus_tags.tsv
```

Existing gene IDs retain their exact tags. Tags for genes absent from the input
remain reserved in the new mapping. New tags are allocated above the highest
reserved or active numeric suffix; alphanumeric suffixes are preserved but do
not contribute a numeric high-water mark. Allocation does not restart by sequence
or gene type. Existing tags are not renumbered to match new coordinate order.
First-time numbering follows the converter's deterministic sequence-ID and
coordinate ordering, not necessarily the original input row order.

DDBJ permits insertion into unused numeric gaps. The tool deliberately does not
automatically fill gaps, because a missing tag may have been retired. You can
explicitly supply an independently verified unused assignment in the map; keep
all historical reservations. Withholding the historical map removes the tool's
ability to protect those historical tags. Gene matching uses stable IDs, not
coordinates or sequence similarity: changed gene IDs require a curated mapping.

Prepare the whole assembly annotation together. Batch mode rejects repeated
prefixes because independent allocations could collide; do not run separate
chromosome jobs with the same prefix and independent maps.

## Feature scope and validation

The gene and its supported transcript/RNA, CDS, exon, intron and UTR descendants
receive the same tag. `repeat_region` never receives an automatically added tag;
a retained input repeat_region already carrying one is rejected for explicit correction.
Other specialist features retain the existing gene-hierarchy policy. Recognising
an RNA component for locus-tag propagation does not guarantee that its feature
term is accepted by the intended Webin build. Orphan RNA features
still require a gene hierarchy or an appropriate source adapter.

One tag cannot be assigned to different gene IDs. Existing input attributes and
the supplied mapping must agree. The preparation report records allocation
policy under `submission_metadata.locus_tag_policy`. Validate the resulting file
with your intended Webin build; local preparation does not establish archive
acceptance.

Reference documentation:

- https://www.ddbj.nig.ac.jp/ddbj/locus_tag-e.html
- https://ena-docs.readthedocs.io/en/latest/submit/study/interactive.html

## Organelle filtering

With `--remove-organelles`, excluded genes receive no new tags. Tags from the
complete supplied mapping remain reserved, including genes now excluded. Existing
tags on excluded features are also reserved when their gene ancestry can be
resolved safely. Keep the cumulative map across filtered and unfiltered builds.
See [ORGANELLES.md](ORGANELLES.md) for matching and exclusion policy.
