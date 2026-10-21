#!/usr/bin/env python3
"""Run the gff3tools validator bundled inside a supplied Webin JAR.

Extract dependencies temporarily and capture Java output in a new log. This is
local, report-only validation: it does not authenticate, upload, check archive
sequence availability or replace the complete Webin submission workflow."""
import argparse,json,subprocess,tempfile,zipfile,re,collections,hashlib
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--jar',required=True);p.add_argument('--gff',required=True);p.add_argument('--log',required=True)
a=p.parse_args()
if Path(a.log).exists():p.error('Choose a new log path')
with tempfile.TemporaryDirectory(prefix='ena-local-validator-') as tmp:
 with zipfile.ZipFile(a.jar) as z:
  for name in z.namelist():
   if name.startswith('BOOT-INF/lib/') and name.endswith('.jar'):
    (Path(tmp)/Path(name).name).write_bytes(z.read(name))
 cmd=['java','-Xmx4g','-cp',str(Path(tmp)/'*'),'uk.ac.ebi.embl.gff3tools.cli.Main','validation']
 cmd+=[str(Path(a.gff).resolve())]
 with open(a.log,'x') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
s=Path(a.log).read_text()
print(json.dumps({'exit_code':r.returncode,'reported_rules':dict(collections.Counter(re.findall(r'Violation of rule (CDS_LENGTH|CDS_INTRON_LENGTH|EXON_LENGTH)\b',s))),'passed_message':bool(re.search(r'The file passed validations',s)),'scope':'Bundled local gff3tools; not authenticated Webin validation or submission'},indent=2))
raise SystemExit(r.returncode)
