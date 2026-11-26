# Contributing

Report reproducible problems through GitHub Issues. Include the toolkit version
or commit, Python version, predictor and output format, command/configuration,
and the relevant error. Use a minimal synthetic or publicly shareable example;
do not include credentials, private paths or unpublished genomic data.

For changes, create a branch, keep modifications focused and add regression tests
for new behavior or bug fixes. Preserve retained prediction coordinates, strands,
CDS phases and product relationships unless a separately documented change is
explicitly intended. Do not invent biological evidence to satisfy a validator.

Run from the repository root:

```bash
python3 -m unittest discover -s tests -v
python3 ena_gffprep.py --help
git diff --check
```

Update relevant documentation and describe the validation performed in the pull
request. Distinguish synthetic tests, real-data preparation and Webin acceptance.
See [code layout](README.md#code-layout) and [compatibility](COMPATIBILITY.md).
