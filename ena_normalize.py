"""Source-independent normalization of explicit feature graphs.

This module never inspects column 2, predictor names or software versions.
Parsing and provenance are separate concerns. Every transformation preserves
coordinates, strand, scores and CDS phases, and refuses ambiguous product merges.
"""
import re
from collections import Counter, defaultdict

from ena_vocabularies import GENE_FEATURES
GENES = GENE_FEATURES
RNAS = {'mRNA','transcript','primary_transcript','rRNA','tRNA','ncRNA',
        'lnc_RNA','miRNA','snRNA','snoRNA','misc_RNA'}
PRODUCT_KEYS = ('protein_id','product','transl_table','exception','pseudo',
                'pseudogene','inference','experiment','partial','locus_tag')
SEGMENT_KEYS = {'ID','original_id','cds_type'}


def _copy(row):
    """Copy mutable row data while retaining its original input coordinates/line."""
    return dict(row, c=list(row['c']), a={k:list(v) for k,v in row['a'].items()})


def _retain_id(attrs, previous):
    """Append an original identifier without discarding older provenance."""
    values = attrs.setdefault('original_id', [])
    if previous is not None and previous not in values:
        values.append(previous)


def _split_leaves(rows, scalar, require, counts, track_original=True):
    """Duplicate multi-parent leaves before transcript-level CDS grouping."""
    nodes = defaultdict(list)
    for row in rows:
        if 'ID' in row['a']: nodes[scalar(row['a'],'ID')].append(row)
    referenced = {p for row in rows for key in ('Parent','Derives_from') for p in row['a'].get(key,[])}
    replacements = {}
    used = set(nodes)
    for fid,parts in nodes.items():
        parents = parts[0]['a'].get('Parent',[])
        if len(parents) < 2: continue
        require(len(set(parents)) == len(parents), f'Duplicate Parent on {fid}')
        require(all(row['a'].get('Parent',[]) == parents for row in parts), f'Inconsistent repeated ID {fid}')
        require(fid not in referenced and all(row['c'][2] not in GENES | RNAS for row in parts), f'Multiple-parent internal or referenced feature {fid}: cannot safely duplicate it')
        require(not any('protein_id' in row['a'] for row in parts), f'Shared feature {fid} carries protein_id; independent product identity needs review')
        for index,parent in enumerate(parents,1):
            new=f'{fid}.ena_parent{index}'
            require(new not in used, 'Generated ID collision: '+new);used.add(new)
        copies=[]
        for row in parts:
            for index,parent in enumerate(parents,1):
                copy=_copy(row);a=copy['a'];a['Parent']=[parent];a['ID']=[f'{fid}.ena_parent{index}']
                if track_original:_retain_id(a,fid)
                copies.append(copy)
        replacements[fid]=copies;counts['shared_features_split']+=1
    result=[];emitted=set()
    for row in rows:
        fid=scalar(row['a'],'ID')
        if fid is None and len(row['a'].get('Parent',[]))>1:
            parents=row['a']['Parent']
            require(len(set(parents))==len(parents), 'Duplicate Parent on anonymous feature')
            for parent in parents:
                copy=_copy(row);copy['a']['Parent']=[parent];result.append(copy)
            counts['shared_features_split']+=1
        elif fid not in replacements:result.append(row)
        elif fid not in emitted:
            result.extend(replacements[fid]);emitted.add(fid)
    return result


def _check_identifiers(rows, scalar, require):
    """Check redundant source identifiers against all reachable parent branches."""
    nodes=defaultdict(list)
    for row in rows:
        if 'ID' in row['a']:nodes[scalar(row['a'],'ID')].append(row)
    for row in rows:
        a=row['a'];fid=scalar(a,'ID');kind=row['c'][2]
        if 'transcript_id' in a and kind not in RNAS:
            pending=list(a.get('Parent',[]));seen=set();rna_ancestors=[]
            while pending:
                parent=pending.pop()
                if parent in seen or parent not in nodes:continue
                seen.add(parent);ancestor=nodes[parent][0]
                if ancestor['c'][2] in RNAS:rna_ancestors.append(ancestor)
                else:pending.extend(ancestor['a'].get('Parent',[]))
            require(len(rna_ancestors)==1 and scalar(a,'transcript_id') in
                    {scalar(rna_ancestors[0]['a'],'ID'), scalar(rna_ancestors[0]['a'],'transcript_id')},
                    f'Line {row["line"]}: transcript_id conflicts with ID/Parent')
        if 'gene_id' not in a:continue
        pending=[row];seen=set();owners=set()
        while pending:
            current=pending.pop()
            if current['c'][2] in GENES:
                owners.add((scalar(current['a'],'ID'), scalar(current['a'],'gene_id')));continue
            for parent in current['a'].get('Parent',[]):
                if parent in seen:continue
                seen.add(parent)
                if parent in nodes:pending.append(nodes[parent][0])
        if owners:
            require(len(owners)==1 and scalar(a,'gene_id') in set(next(iter(owners))),
                    f'Line {row["line"]}: gene_id conflicts with ancestor gene')


def normalize(rows, scalar, require, mode='standard', cds_grouping='auto'):
    """Return normalized rows and rule counts using structure alone.

    Rows are mutated in place; the returned list can also contain copied leaves.
    auto merges compatible uniquely identified segments, preserves already-grouped
    features and distinct explicitly identified proteins, and refuses ambiguity.
    single-product explicitly requests one compatible CDS per transcript; it never
    overrides conflicting protein/evidence metadata. preserve keeps CDS group IDs.
    strict skips type/phase/CDS rewrites; single-parent leaf expansion and final
    graph checks still apply.
    """
    require(mode in {'standard','strict'}, 'Invalid normalization mode')
    require(cds_grouping in {'auto','single-product','preserve'}, 'Invalid CDS grouping policy')
    counts=Counter()
    rows=_split_leaves(rows,scalar,require,counts,track_original=mode=='standard')
    if mode=='strict':return rows,dict(counts)
    _check_identifiers(rows,scalar,require)
    coding={p for row in rows if row['c'][2]=='CDS' for p in row['a'].get('Parent',[])}
    groups=defaultdict(list)
    used={scalar(row['a'],'ID') for row in rows if 'ID' in row['a']}
    for row in rows:
        c,a=row['c'],row['a'];fid=scalar(a,'ID')
        if c[2]!='CDS' and c[7] in {'0','1','2'}:
            require('original_phase' not in a or a['original_phase']==[c[7]], f'Line {row["line"]}: original_phase conflict')
            a['original_phase']=[c[7]];c[7]='.';row['normalized_phase']=True;counts['non_CDS_phases_normalized']+=1
        if c[2]=='transcript' and fid in coding:
            c[2]='mRNA';counts['coding_transcripts_to_mRNA']+=1
        if 'transcript_id' in a and c[2] not in RNAS:
            require('original_transcript_id' not in a or a['original_transcript_id']==a['transcript_id'], 'Conflicting original_transcript_id')
            a['original_transcript_id']=a.pop('transcript_id');counts['transcript_id_attributes_mapped']+=1
        if c[2]=='CDS':
            parents=a.get('Parent',[])
            require(len(parents)==1, f'Line {row["line"]}: Orphan or shared anonymous CDS requires one transcript Parent')
            groups[parents[0]].append(row)
    transcripts={scalar(row['a'],'ID') for row in rows if row['c'][2]=='mRNA'}
    all_ids={scalar(row['a'],'ID') for row in rows if 'ID' in row['a']}
    refs={p for row in rows for key in ('Parent','Derives_from') for p in row['a'].get(key,[])}
    for parent,parts in groups.items():
        require(parent in all_ids, f'Missing Parent {parent} for CDS')
        counts['cds_transcripts']+=1
        old={scalar(row['a'],'ID') for row in parts}
        signatures={(row['c'][0],row['c'][6]) for row in parts}
        require(len(signatures)==1, f'CDS sequence/strand mismatch under {parent}')
        # One valid existing logical group remains unchanged. The numbered
        # singleton convention is normalized structurally for legacy compatibility.
        numbered=parent in transcripts and all(fid is not None and re.fullmatch(re.escape(parent)+r'\.CDS\.[0-9]+',fid) for fid in old)
        if cds_grouping=='preserve':
            require(None not in old, f'CDS grouping preserve requires explicit IDs: {parent}')
            counts['cds_groups_preserved']+=len(old);continue
        if len(old)==1 and None not in old and not numbered:
            counts['cds_groups_preserved']+=1;continue
        by_id=defaultdict(list)
        for row in parts:by_id[scalar(row['a'],'ID')].append(row)
        proteins=[{scalar(row['a'],'protein_id') for row in group} for group in by_id.values()]
        distinct_proteins=(None not in old and len(old)>1 and all(len(values)==1 and None not in values for values in proteins)
                           and len({next(iter(values)) for values in proteins})==len(proteins))
        if cds_grouping=='auto' and distinct_proteins:
            counts['explicit_product_groups_preserved']+=len(old);continue
        require(parent in transcripts, f'CDS regrouping requires one coding transcript parent: {parent}; preserve existing groups or correct the hierarchy')
        require(None not in old or old=={None}, f'Mixed identified/anonymous CDS segments require explicit correction: {parent}')
        require(cds_grouping=='single-product' or old=={None} or len(old)==len(parts),
                f'Ambiguous existing multipart CDS groups for {parent}; use preserve or explicitly confirm single-product')
        require(not (old-{None}) & refs, f'Referenced CDS segment IDs require explicit review: {parent}')
        intervals=sorted((row['start'],row['end']) for row in parts)
        require(all(a[1]<b[0] for a,b in zip(intervals,intervals[1:])), f'Overlapping CDS products for {parent}')
        for key in PRODUCT_KEYS:
            values={tuple(row['a'].get(key,[])) for row in parts}
            require(len(values)==1, f'Conflicting {key} in CDS attributes across segments: {parent}')
        metadata=[{key:value for key,value in row['a'].items() if key not in SEGMENT_KEYS} for row in parts]
        require(all(a==metadata[0] for a in metadata), f'Conflicting CDS attributes across segments: {parent}; per-segment metadata needs explicit review')
        new=parent+'.CDS'
        require(new not in used, 'Normalized CDS ID collision: '+new);used.add(new)
        for row in parts:
            previous=scalar(row['a'],'ID')
            if previous is not None:_retain_id(row['a'],previous)
            row['a']['ID']=[new];counts['cds_rows_renamed' if previous else 'cds_rows_assigned_ID']+=1
        counts['cds_groups_normalized']+=1
    return rows,dict(counts)
