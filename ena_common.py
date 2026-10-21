"""Shared CLI validation, text I/O and package-publication helpers.

These helpers create local artifacts only. They do not contact ENA, manage
credentials or establish biological validity."""
import csv
import gzip
import json
import re
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path

VERSION = '0.5.5'

def check(condition, message):
    """Raise ValueError with message when condition is false."""
    if not condition:
        raise ValueError(message)

def text_field(value, name):
    """Return a nonempty printable ASCII metadata value without tabs or newlines."""
    check(isinstance(value, str) and value.strip() and all(32 <= ord(c) < 127 for c in value),
          f'{name} must be nonempty printable ASCII, without tabs/newlines')
    return value

def safe_label(value):
    """Validate a package label suitable for a single filename component."""
    check(re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', value or '') is not None,
          'Label must start with a letter/digit and contain only letters, digits, dots, underscores or hyphens')
    return value

def input_stem(path):
    """Derive a sanitised default label by removing recognised input suffixes."""
    value = Path(path).name
    for suffix in ['.gz', '.gff3', '.gff', '.gtf', '.fasta', '.fna', '.fa']:
        if value.lower().endswith(suffix):
            value = value[:-len(suffix)]
    return re.sub(r'[^A-Za-z0-9._-]+', '_', value).strip('._-') or 'annotation'

def open_text(path, mode='rt'):
    """Open plain or .gz-suffixed text using ASCII encoding."""
    return gzip.open(path, mode, encoding='ascii') if str(path).endswith('.gz') else open(path, mode, encoding='ascii')

def tsv_rows(path, allowed, required):
    """Read strict TSV rows with allowed and required column names.
    
    Strip cell whitespace, skip blank rows and reject malformed row widths.
    Return dictionaries; blank cells remain empty strings for caller-side inheritance."""
    with open(path, newline='', encoding='utf-8-sig') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        check(len(fields) == len(set(fields)), 'Duplicate TSV column names')
        check(set(required) <= set(fields), 'Missing TSV columns: ' + ', '.join(sorted(set(required)-set(fields))))
        check(not set(fields)-set(allowed), 'Unknown TSV columns: ' + ', '.join(sorted(set(fields)-set(allowed))))
        result = []
        for number, row in enumerate(reader, 2):
            check(None not in row and all(v is not None for v in row.values()), f'TSV row {number}: wrong number of columns')
            if any(v.strip() for v in row.values()):
                result.append({k:v.strip() for k,v in row.items()})
        return result

@contextmanager
def new_package(destination):
    """Yield a private staging directory, then publish it as a new destination.
    
    Refuse an existing destination and remove staging files if preparation fails.
    The staging directory shares the destination filesystem so rename publishes
    the completed package only after its writers have closed."""
    target = Path(destination).resolve()
    check(not target.exists(), f'Output directory already exists: {target}')
    target.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.ena-prep-', dir=target.parent))
    try:
        yield staging
        check(not target.exists(), f'Output directory appeared during preparation: {target}')
        # rename is a same-filesystem publication; all files are closed before this point.
        staging.rename(target)
    finally:
        if staging.exists():
            shutil.rmtree(staging)

def write_json(path, value):
    """Write indented JSON with ASCII escaping and a trailing newline."""
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=True)+'\n', encoding='ascii')

def write_manifest(path, fields):
    """Write tab-separated manifest fields after validating their text values."""
    with open(path, 'w', encoding='ascii') as handle:
        for key, value in fields:
            handle.write(key+'\t'+text_field(str(value), key)+'\n')

def validation_instructions(directory, context):
    """Write VALIDATE.txt for a prepared package and the given Webin context.
    
    This creates instructions only; it neither handles credentials nor runs Webin."""
    Path(directory, 'VALIDATE.txt').write_text(
        'Change into this output directory before running the supplied pilot JAR.\n'
        'Set WEBIN_USERNAME and WEBIN_PASSWORD in your environment.\n'
        'Follow your ENA contact\'s routing instructions; the supplied examples omit -test.\n\n'
        'java -Xmx4g -XX:ActiveProcessorCount=2 -jar /absolute/path/to/webin-cli-10.0.0-beta.jar \\\n'
        f'  -context {context} -manifest manifest.txt \\\n'
        '  -username "$WEBIN_USERNAME" -passwordEnv WEBIN_PASSWORD -validate\n\n'
        'This command validates only. No toolkit command submits data.\n', encoding='ascii')
