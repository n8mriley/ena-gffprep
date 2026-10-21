#!/usr/bin/env python3
"""Add explicit CDS prediction provenance to a separate GFF3 output.

Use only evidence text that accurately describes the prediction method. An
inference attribute does not prove that a short CDS or other prediction is valid."""
import argparse,gzip
from pathlib import Path
from urllib.parse import quote
from ena_gffprep import staged_outputs

def add(source,destination,value):
 """Write a new GFF3 with the supplied inference on CDS rows; return their count.
 
 Require explicit printable evidence text and refuse existing CDS inference
 rather than overwriting it. Escape the attribute value and stage the output.
 This helper records provenance; it does not establish biological correctness."""
 if Path(source).resolve()==Path(destination).resolve():raise ValueError('Choose a new output')
 if Path(destination).exists():raise ValueError('Output already exists')
 if not value or any(ord(c)<32 or ord(c)>126 for c in value):raise ValueError('Inference must be nonempty ASCII without controls')
 op=gzip.open if str(source).endswith('.gz') else open
 out=gzip.open if str(destination).endswith('.gz') else open
 changed=0
 with staged_outputs([destination]) as staged, op(source,'rt',encoding='ascii') as f:
  with out(staged[destination],'wt',encoding='ascii') as g:
   for l in f:
    if l.strip() and not l.startswith('#'):
     c=l.rstrip('\r\n').split('\t')
     if len(c)!=9:raise ValueError('Expected nine columns')
     if c[2]=='CDS':
      if any(x.startswith('inference=') for x in c[8].split(';')):raise ValueError('Existing inference: use a dedicated merge review')
      c[8]=('' if c[8]=='.' else c[8].rstrip(';')+';')+'inference='+quote(value,safe=':._-+|/@()[] ')
      l='\t'.join(c)+'\n';changed+=1
    g.write(l)
 return changed
if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__)
 for k in ['input','output','inference']:p.add_argument('--'+k,required=True)
 a=p.parse_args();print('CDS rows annotated:',add(a.input,a.output,a.inference))
