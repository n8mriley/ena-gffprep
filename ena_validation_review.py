#!/usr/bin/env python3
"""Associate supported local validator diagnostics with prepared GFF3 features.

Keep the exact validated GFF3 alongside its log because associations rely on
physical line numbers. Additional coordinate-based QC is advisory and does not
translate sequences or alter predictions."""
import argparse, collections, csv, gzip, hashlib, json, re
from pathlib import Path
from urllib.parse import unquote

def review(gff, log):
    """Return summary, validator diagnostics, short-CDS and short-intron tables.
    
    Use the exact GFF3 validated by the log: physical line numbers associate
    diagnostics with features. Recognise the supported violation-log syntax only.
    Coordinate-only advisory thresholds do not demonstrate biological errors."""
    rows, nodes, cds = {}, {}, collections.defaultdict(list)
    op = gzip.open if str(gff).endswith('.gz') else open
    with op(gff, 'rt', encoding='ascii') as f:
        for n, line in enumerate(f, 1):
            if not line.strip() or line.startswith('#'): continue
            c = line.rstrip('\r\n').split('\t')
            if len(c) != 9: raise ValueError('Expected nine columns at line %s' % n)
            a = {k:[unquote(v) for v in vs.split(',')] for k,vs in (x.split('=',1) for x in c[8].rstrip(';').split(';') if '=' in x)}
            r = dict(line=n, seqid=unquote(c[0]), type=c[2], start=int(c[3]), end=int(c[4]), strand=c[6], phase=c[7], attrs=a)
            rows[n] = r
            if 'ID' in a: nodes.setdefault(a['ID'][0], r)
            if c[2] == 'CDS':
                if 'ID' not in a: raise ValueError('Use prepared GFF3 with logical CDS IDs')
                cds[a['ID'][0]].append(r)
    def identity(r):
        """Resolve row identifiers and ancestor gene metadata through single parents.
        
        Detect cycles and preserve unresolved metadata rather than inventing matches."""
        a = r['attrs']; parent = a.get('Parent',[''])[0]; gene = ''; cursor = r; seen=set()
        while cursor:
            ca=cursor['attrs']; fid=ca.get('ID',[''])[0]
            if fid in seen: raise ValueError('Parent cycle')
            seen.add(fid)
            if cursor['type']=='gene': gene=fid; break
            cursor=nodes.get(ca.get('Parent',[''])[0])
        return dict(seqid=r['seqid'], strand=r['strand'], feature_id=a.get('ID',[''])[0], parent_id=parent, gene_id=gene, locus_tag=a.get('locus_tag',[''])[0])
    diagnostics=[]
    pattern=re.compile(r'^(ERROR|WARNING): Violation of rule (\w+).*?on line (\d+)(.*)$')
    for line in Path(log).read_text().splitlines():
        m=pattern.match(line)
        if not m: continue
        n=int(m[3])
        if n not in rows: raise ValueError('Diagnostic line %s absent from GFF3; confirm exact validated file' % n)
        r=rows[n]
        accession=re.search(r'accession "([^"]+)"',line)
        if accession and accession[1]!=r['seqid']: raise ValueError('Diagnostic accession mismatch at line %s' % n)
        diagnostics.append(dict(severity=m[1], rule=m[2], gff_line=n, **identity(r), start=r['start'], end=r['end'], interval_nt=r['end']-r['start']+1, message=line))
    short=[]; introns=[]
    reported_short={d['feature_id'] for d in diagnostics if d['rule']=='CDS_LENGTH'}
    for fid,parts in cds.items():
        p=sorted(parts,key=lambda r:r['start']); first=p[-1] if p[0]['strand']=='-' else p[0]
        length=sum(r['end']-r['start']+1 for r in p)-int(first['phase'])
        if length<90:
            short.append(dict(**identity(first), first_gff_line=min(r['line'] for r in p), coding_nt=length, segments=len(p), reported_CDS_LENGTH=fid in reported_short, candidate_at_or_below_75_nt=length<=75, partial=','.join(sorted({v for r in p for v in r['attrs'].get('partial',[])})), has_evidence_or_inference=any('experiment' in r['attrs'] or 'inference' in r['attrs'] for r in p), has_exception_or_pseudo=any(any(k in r['attrs'] for k in ['exception','pseudo','pseudogene']) for r in p)))
        for left,right in zip(p,p[1:]):
            gap=right['start']-left['end']-1
            if gap<10:
                introns.append(dict(**identity(first), left_gff_line=left['line'], right_gff_line=right['line'], intron_start=left['end']+1, intron_end=right['start']-1, intron_nt=gap))
    summary=dict(diagnostic_counts=dict(collections.Counter(d['severity']+':'+d['rule'] for d in diagnostics)), logical_CDS=len(cds), CDS_under_90_nt=len(short), CDS_at_or_below_75_nt=sum(s['coding_nt']<=75 for s in short), CDS_introns_under_10_nt=len(introns), CDSs_with_introns_under_10_nt=len({r['feature_id'] for r in introns}), source_sha256={str(p):hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in [gff,log]}, limitations=['Report line mapping requires the exact validated GFF3; accession checks detect only some mismatches.', '75 nt is a review boundary inferred from observed diagnostics, not a complete reproduction of validator translation rules.', '90 nt remains a descriptive review threshold.', 'No FASTA: translation, terminal stops, splice motifs and biological completeness are unverified.', 'Intron calculation assumes ordinary same-sequence collinear non-overlapping CDS intervals; specialised biology needs separate review.', 'No annotation changes or submission performed.'])
    return summary, diagnostics, short, introns

def main():
    """Write validator-review TSV and JSON reports into a new output directory.
    
    The directory is not staged as a transaction; failures may leave partial output."""
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--gff',required=True);p.add_argument('--webin-report',required=True);p.add_argument('--outdir',required=True);a=p.parse_args()
    summary,diagnostics,short,introns=review(a.gff,a.webin_report)
    out=Path(a.outdir);out.mkdir(parents=True,exist_ok=False)
    for name,data in [('webin_findings',diagnostics),('short_cds_review',short),('short_cds_introns',introns)]:
        with (out/(name+'.tsv')).open('w',newline='') as f:
            if data:
                w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
