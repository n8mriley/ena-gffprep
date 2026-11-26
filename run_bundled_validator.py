#!/usr/bin/env python3
"""Run bundled gff3tools locally without fixes, authentication or submission."""
import argparse
import collections
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile
import zipfile


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('jar', 'gff', 'log'):
        p.add_argument('--'+key, required=True)
    p.add_argument('--sequence')
    p.add_argument('--metrics')
    p.add_argument('--json-report')
    p.add_argument('--memory', default='4g')
    p.add_argument('--threads', type=int, default=2)
    a = p.parse_args()
    if a.threads < 1 or not re.fullmatch(r'[1-9][0-9]*[kKmMgG]?', a.memory):
        p.error('Invalid JVM resources')
    outputs = [Path(x).resolve() for x in (a.log, a.metrics, a.json_report) if x]
    if len(set(outputs)) != len(outputs) or any(x.exists() for x in outputs):
        p.error('Choose distinct new output paths')
    with tempfile.TemporaryDirectory(prefix='ena-local-validator-') as tmp:
        with zipfile.ZipFile(a.jar) as z:
            libs = [n for n in z.namelist() if n.startswith('BOOT-INF/lib/') and n.endswith('.jar')]
            validators = [Path(n).name for n in libs if Path(n).name.startswith('gff3tools-')]
            if len(validators) != 1:
                p.error('Expected exactly one bundled gff3tools library')
            for name in libs:
                (Path(tmp)/Path(name).name).write_bytes(z.read(name))
        cmd = ['java', '-Xmx'+a.memory, '-XX:ActiveProcessorCount='+str(a.threads), '-cp', str(Path(tmp)/'*'), 'uk.ac.ebi.embl.gff3tools.cli.Main', 'validation']
        if a.sequence:
            cmd += ['--sequence', str(Path(a.sequence).resolve())]
        if a.metrics:
            cmd += ['--metrics', str(Path(a.metrics).resolve()), '--metrics-format', 'json']
        cmd += [str(Path(a.gff).resolve())]
        with open(a.log, 'x') as log:
            result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    text = Path(a.log).read_text(errors='replace')
    report = {'exit_code': result.returncode, 'reported_rules': dict(collections.Counter(re.findall(r'Violation of rule ([A-Z][A-Z0-9_]+)\b', text))), 'passed_message': bool(re.search(r'The file (?:has )?passed (?:all )?validations', text)), 'scope': 'Bundled standalone gff3tools; not authenticated Webin validation or submission', 'validator': validators[0], 'jar_sha256': hashlib.sha256(Path(a.jar).read_bytes()).hexdigest(), 'java_version': subprocess.run(['java','-version'], capture_output=True, text=True).stderr.strip(), 'command': cmd, 'sequence': a.sequence, 'metrics': a.metrics}
    rendered = json.dumps(report, indent=2)
    if a.json_report:
        with open(a.json_report, 'x') as f:
            f.write(rendered+'\n')
    print(rendered)
    return result.returncode


if __name__ == '__main__':
    raise SystemExit(main())
