# Examples

The accession numbers, prefixes and annotation paths in this directory are
placeholders. `pilot.gff3` is synthetic. Do not submit these examples.

| File | Use |
| --- | --- |
| `pilot.gff3` | Small structural fixture with known local bounds and no real prediction provenance. |
| `nuclear_defaults.json` | Shared translation-table fallback and output root; paths are relative to the JSON file. |
| `annotations.tsv` | Multi-annotator batch template; replace paths, accessions and prefixes before use. |

## Run the synthetic fixture

From the toolkit directory, choose a new output directory:

```bash
python3 ena_gffprep.py \
  --input examples/pilot.gff3 \
  --assembly-accession GCA_123456789.1 \
  --bioproject PRJEB123456 \
  --locus-tag-prefix EXAMPLE \
  --no-inference --no-ena-lengths \
  --output-dir example-preparation
```

The explicit no-inference policy avoids assigning invented provenance to this
fixture. Its sequence-region declaration supplies bounds, so no lookup is needed.
The accessions/prefix are illustrative, not verified registration or archive links.

## Adapt the batch template

The table includes direct-source and experimental evidence-based pipeline examples. BRAKER compatibility is not fully tested; see the [scope statement](../README.md#scope-and-maturity). Replace each `ACTUAL_VERSION`
with the version used. Empty cells inherit config/CLI defaults; they do not clear
a default. Keep predictor/model defaults out of shared settings when jobs differ.
Each job needs its own registered prefix or a deliberately coordinated allocation.

```bash
python3 ena_gffprep.py \
  --config examples/nuclear_defaults.json \
  --batch examples/annotations.tsv \
  --batch-summary example-batch-summary.json
```

Paths in row cells are relative to the TSV unless absolute. Add optional columns
such as `fai`, `locus_tag_map`, `remove_organelles` or `sequence_metadata` as needed.
TSV booleans use `true`/`false`. See the [main guide](../README.md#configuration-and-batches)
for precedence and batch restrictions.
