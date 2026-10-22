# Upgrade an existing Git checkout to 0.6.0

Extract this distribution outside your existing repository. Keep the complete
script directory, including the new ena_organelles.py module. The release also
includes the current Annotator / Version / Model description format and custom
description behaviour. Existing manifests are not edited automatically.

From your repository, check for local changes and create an update branch:

```bash
git status
git switch main
git pull --ff-only
git switch -c add-organelle-filter
```

Preserve any uncommitted work before copying. Replace the source path below with
your extracted package location. Preview, then copy:

```bash
rsync -avnc --exclude='SHA256SUMS' --exclude='__pycache__/' \
  /path/to/ena_gffprep-0.6.0/ ./
rsync -avc --exclude='SHA256SUMS' --exclude='__pycache__/' \
  /path/to/ena_gffprep-0.6.0/ ./
```

Do not add --delete. These commands retain your .git history, .github workflow,
.gitignore, licence and local Webin JAR. Package files such as README.md are
updated; review custom local documentation after copying.

```bash
git diff --check
git diff --stat
python3 ena_gffprep.py --version
python3 -m unittest discover -s tests -q
```

Expect version 0.6.0 and 120 passing tests. Stage intended files explicitly,
leaving local backups and real submission data untracked:

```bash
git add ena_gffprep.py ena_gffprep_core.py ena_common.py ena_organelles.py \
  tests/test_organelles.py tests/test_predictor_model.py \
  README.md ORGANELLES.md UPGRADING.md CHANGELOG.md VALIDATION.md release_verification.json
git diff --cached --stat
git commit -m "Add optional NCBI-based organelle annotation exclusion"
git push -u origin add-organelle-filter
```

Create a pull request into main, check the automated test result and merge. Then:

```bash
git switch main
git pull --ff-only
```

Use a new output directory when trying --remove-organelles. For annotation
updates supply the complete previous locus-tag map; excluded genes stay reserved.
See ORGANELLES.md for online, cached and offline usage.
