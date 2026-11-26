# Local reference data

The v0.9.0 development branch requires two locally generated indexes for
annotation preparation and the full regression suite. Neither index is shipped.
Python 3.9+ is sufficient for setup; Java is only needed to execute a validator.
Obtain external materials under their applicable terms (see
[third-party notices](THIRD_PARTY_NOTICES.md)).

## 1. Download Sequence Ontology

Run from the repository root. Keep downloaded resources outside the checkout.

```bash
mkdir -p ../ena-reference-data
wget -O ../ena-reference-data/so.obo \
  https://raw.githubusercontent.com/The-Sequence-Ontology/SO-Ontologies/4340c143bac3578bba0c013d02e8c5f0e51ad14a/Ontology_Files/so.obo
sha256sum ../ena-reference-data/so.obo
```

Expected SHA-256: `22a8f3ec2b49125dbb6cee8456f0bca86dd8c98f433165ffa4b554da4f155204`.
The pinned revision contains data version `2026-08-07`. Do not substitute the
moving latest release or change the checksum to bypass a mismatch.

## 2. Obtain compatible ENA feature mappings

Download the Java JAR from the [official Webin CLI releases](https://github.com/enasequence/webin-cli/releases).
The public Webin CLI 11.0.0 release was verified to contain the pinned embedded gff3tools bytes.
Save the downloaded JAR outside the repository and set its actual path:

```bash
wget -O ../ena-reference-data/webin-cli-11.0.0.jar \
  https://github.com/enasequence/webin-cli/releases/download/11.0.0/webin-cli-11.0.0.jar

export WEBIN_JAR="$PWD/../ena-reference-data/webin-cli-11.0.0.jar"
```

Setup accepts either a Webin CLI JAR containing exactly one embedded gff3tools
JAR or the standalone gff3tools JAR itself. The embedded/standalone gff3tools
5.5.1 bytes must hash to
`722b23c530447623597b58f5d3a929bcffe953062480ff3d2fddd6fc864c66cc`.
The outer Webin JAR checksum is not the setup pin. Check a candidate without
running Java or extracting files to disk:

```bash
python3 - "$WEBIN_JAR" <<'PY'
import hashlib, sys, zipfile
from pathlib import Path
p = Path(sys.argv[1])
with zipfile.ZipFile(p) as z:
    if 'feature-mapping.tsv' in z.namelist():
        raw = p.read_bytes()
    else:
        names = [n for n in z.namelist()
                 if n.startswith('BOOT-INF/lib/gff3tools-') and n.endswith('.jar')]
        if len(names) != 1:
            raise SystemExit('Expected exactly one embedded gff3tools JAR')
        raw = z.read(names[0])
print(hashlib.sha256(raw).hexdigest())
PY
```

If the hash differs, obtain the compatible build from ENA or contact the project
maintainers. Do not edit the pin. The public release was verified on the development HPC. A newer Webin build may still
be appropriate for final submission validation (see [Validation](VALIDATION.md)).

## 3. Acquire the Dbxref registry

Open the [official INSDC Dbxref vocabulary](https://www.insdc.org/submitting-standards/dbxref-qualifier-vocabulary/).
Download the source page for automated extraction and provenance:

```bash
wget -O ../ena-reference-data/dbxref-registry.html \
  https://www.insdc.org/submitting-standards/dbxref-qualifier-vocabulary/
```

Extract the example tags and review the reported discrepancies:

```bash
python3 tools/extract_dbxref_tags.py ../ena-reference-data/dbxref-registry.html ../ena-reference-data/dbxref-tags.txt
# After reviewing discrepancies, rerun with explicit acknowledgement:
python3 tools/extract_dbxref_tags.py ../ena-reference-data/dbxref-registry.html ../ena-reference-data/dbxref-tags.txt --accept-example-tags
```

The first command refuses to write if discrepancies exist; neither command overwrites
an existing list. For the HTML retrieved on 2026-11-26, the source SHA-256 was
`5707672eb39c2333eb950e28c07c4905495573c3ad8c26a2e31fa1eccd34e1fc` and the
126-tag extracted list SHA-256 was
`e57100ead84f0f2b1d8c070f4582ce2725365d17b3feafad8d089d85fd18d4d2`.
The preparation policy accepts both BOLD/Bold and InterPro/Interpro when either
variant is present, excludes the apparent typo Phytozme, and retains Phytozome
and TAIR if present. This gives 127 effective tags for that snapshot. No input
attributes are rewritten. The list remains unmodified and case-sensitive.
A changed upstream registry requires review; setup checks format, not completeness.

## 4. Generate and check the indexes

```bash
python3 tools/setup_references.py \
  --obo ../ena-reference-data/so.obo \
  --jar "$WEBIN_JAR" \
  --dbxref-tags ../ena-reference-data/dbxref-tags.txt

# In a Git checkout, verify that generated indexes are ignored:
git check-ignore -v so_snapshot.json local_vocabularies.json
python3 -m unittest discover -s tests -v
```

Expected setup message: `Local reference indexes generated; do not commit them.`
Both JSON files are written beside the Python modules. The full suite passed 259
tests in an independent clean installation using freshly downloaded public
references on 2026-11-26. Setup
itself does not access the network. Repeat it for each new checkout/worktree.
Keep the same inputs for repeatable preparation; retain their checksums privately.
The preparation JSON records the ontology version and source hash, gff3tools
hash and Dbxref list hash. The list is not checksum-pinned across installations.

Missing indexes require setup, not disabling validation. A checksum mismatch
requires the correct upstream file. An unknown Dbxref tag requires checking the
registry and regenerating the index if appropriate, not inventing a tag.
Do not commit or publish downloaded resources or generated indexes. Network-free
preparation additionally requires local sequence lengths/metadata or disabled
ENA lookup where appropriate; reference setup alone does not supply assembly metadata.
