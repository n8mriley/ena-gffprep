#!/usr/bin/env python3
"""Prepare reproducible ENA decoupled-annotation packages.

The CLI resolves configuration, predictor provenance and output naming before
delegating GFF3 representation changes to ena_gffprep_core. Predictions are
preserved; preparation is separate from official Webin validation and submission."""
import argparse
import csv
import json
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import ena_gffprep_core as core
# Preserve the programmatic preparation API used by earlier callers and QC helpers.
from ena_gffprep_core import *
from ena_common import (VERSION, check, input_stem, new_package, open_text, safe_label,
                        text_field, tsv_rows, validation_instructions, write_json)

PATH_KEYS = {'input', 'fai', 'assembly_report', 'locus_tag_map', 'output_dir', 'out_root',
             'output', 'report', 'manifest', 'write_locus_tag_map'}
SETTING_KEYS = PATH_KEYS | {'source_format', 'assembly_accession', 'taxon_id', 'locus_tag_prefix',
    'translation_table', 'predictor', 'predictor_version', 'predictor_model', 'inference', 'label', 'name',
    'description', 'accession_column', 'width', 'start', 'step', 'default_product',
    'no_inference', 'strict_warnings', 'use_annevo_lengths'}


def detect_source(path):
    """Return (predictor, version) from recognised leading comments.
    
    ANNEVO headers may identify the predictor without specifying a version.
    Only ANNEVO is detected automatically. Conflicting versions raise ValueError;
    filenames are never treated as evidence of prediction provenance."""
    versions = set()
    recognised = False
    with open_text(path) as handle:
        for number, line in enumerate(handle):
            if line.startswith('#') and re.search(r'\b(?:generated\s+(?:with|by)\s+ANNEVO|ANNEVO\s+is\s+an?\s+ab\s+initio)', line, re.I):
                recognised = True
            match = re.search(r'\bANNEVO\s+(?:version\s*|v\s*)?([0-9]+(?:\.[0-9]+)+)', line, re.I) if line.startswith('#') else None
            if match:
                versions.add(match[1])
            if line.strip() and not line.startswith('#'):
                break
            if number > 10000:
                break
    check(len(versions) <= 1, 'Conflicting ANNEVO versions in input header')
    return ('ANNEVO', next(iter(versions), None)) if versions or recognised else (None, None)


def configuration(path):
    """Load validated JSON defaults, resolving file paths relative to the config.
    
    Return an empty dictionary when no configuration file is supplied."""
    if not path:
        return {}
    data = json.loads(Path(path).read_text())
    check(isinstance(data, dict) and not set(data)-SETTING_KEYS, 'Config must be an object containing documented settings only')
    return relative_paths(data, Path(path).resolve().parent)


def relative_paths(data, base):
    """Return a copy with known relative path settings resolved against base.
    
    Metadata strings, including descriptions and accessions, are left untouched."""
    result = dict(data)
    for key in PATH_KEYS:
        if result.get(key) and not Path(result[key]).is_absolute():
            result[key] = str((base / result[key]).resolve())
    return result


def resolve_settings(settings):
    """Validate and resolve one annotation job without writing its output files.
    
    Derive labels, provenance and conservative defaults. Taxonomy is optional:
    the core preserves a valid input species directive when no taxon is supplied.
    Translation table 1 is a fallback assumption, not an organelle code detector."""
    s = {k:v for k,v in settings.items() if v is not None and v != ''}
    check(s.get('input'), '--input is required')
    check(Path(s['input']).is_file(), f'Input file not found: {s["input"]}')
    detected, version = detect_source(s['input'])
    source_format = s.get('source_format', 'auto')
    check(source_format in {'auto', 'annevo', 'generic'}, 'Unsupported source_format')
    s['source_format'] = ('annevo' if (s.get('predictor') or detected or '').casefold() == 'annevo' else 'generic') if source_format == 'auto' else source_format
    if detected:
        check(not s.get('predictor') or s['predictor'].lower() == detected.lower(), 'Predictor conflicts with ANNEVO header')
        check(version is None or not s.get('predictor_version') or str(s['predictor_version']) == version, 'Predictor version conflicts with input header')
        s.setdefault('predictor', detected)
        if version is not None:
            s.setdefault('predictor_version', version)
    # Inference records the actual prediction method, not experimental support.
    if s.get('no_inference'):
        check(not s.get('inference'), '--no-inference cannot be combined with --inference')
        s['inference'] = None
    elif not s.get('inference'):
        check(s.get('predictor') and s.get('predictor_version'),
              'Cannot determine predictor/version. Supply both --predictor and --predictor-version, --inference, or explicit --no-inference')
        s['inference'] = f'ab initio prediction:{s["predictor"]}:{s["predictor_version"]}'
    elif detected and s.get('predictor_version'):
        check(s['inference'] == f"ab initio prediction:{detected}:{s['predictor_version']}",
              'Explicit inference conflicts with detected ANNEVO provenance')
    # Taxonomy is optional; the accession and registered tag namespace are required.
    for key in ['assembly_accession', 'locus_tag_prefix']:
        check(s.get(key), '--'+key.replace('_','-')+' is required (CLI, config or batch row)')
    check(re.fullmatch(r'GC[AF]_\d+\.\d+', s['assembly_accession']) is not None, 'Use a GCA/GCF assembly accession with version')
    for key, default in [('translation_table',1), ('width',6), ('start',1), ('step',1)]:
        s[key] = int(s.get(key, default))
    s['taxon_id'] = int(s['taxon_id']) if s.get('taxon_id') is not None else None
    s.setdefault('accession_column','genbank')
    s.setdefault('strict_warnings',False)
    s.setdefault('default_product',None)
    s.setdefault('use_annevo_lengths',False)
    s['automatic_annevo_lengths'] = bool(s['source_format']=='annevo' and not s.get('fai') and not s.get('assembly_report'))
    s['use_annevo_lengths'] = s['use_annevo_lengths'] or s['automatic_annevo_lengths']
    if s.get('predictor_model'):
        text_field(s['predictor_model'], 'predictor_model')
        check(s.get('predictor') and s.get('predictor_version'),
              '--predictor-model requires predictor and predictor-version (explicit or detected)')
    if s.get('predictor'):
        # Preserve the requested spelling; normalise unsafe filename characters only.
        components = [s['assembly_accession'], s['predictor']]
        if s.get('predictor_version'):
            components.append(str(s['predictor_version']))
        if s.get('predictor_model'):
            components.append(s['predictor_model'])
        label = '.'.join(re.sub(r'[^A-Za-z0-9._-]+', '_', c) for c in components)
    else:
        label = s['assembly_accession']+'.'+input_stem(s['input'])
    s['label'] = safe_label(s.get('label', label))
    s['output_dir'] = str(Path(s.get('output_dir') or Path(s.get('out_root','prepared_annotations'))/s['label']).resolve())
    s.setdefault('name',s['label']+' decoupled annotation')
    s.setdefault('description','Ab initio genome annotation prepared for ENA decoupled submission.' +
                 (' Prediction provenance: '+s['inference']+'.' if s.get('inference') else ' Existing evidence metadata preserved.'))
    if s.get('predictor_model'):
        # Always retain model provenance, including with a custom description.
        s['description'] += ' Predictor model: '+s['predictor_model']+'.'
    return s


def prepare_package(settings):
    """Create a new annotation package and return a compact outcome dictionary.
    
    Stage the GFF3, locus-tag map, manifest and reports together before publication.
    Existing destination directories are refused. Biological QC remains advisory;
    this function does not authenticate with Webin or submit anything."""
    s = resolve_settings(settings)
    label = s['label']
    with new_package(s['output_dir']) as directory:
        values = dict(s)
        for key in ['fai','assembly_report','locus_tag_map']:
            values.setdefault(key,None)
        values.update(output=str(directory/(label+'.ena.gff3.gz')),
                      write_locus_tag_map=str(directory/'locus_tags.tsv'),
                      report=str(directory/'preparation.json'),manifest=str(directory/'manifest.txt'))
        report = core.prepare(SimpleNamespace(**values))
        # Keep core's hashes under stable source names, never transient output paths.
        report['workflow'] = {'toolkit_version':VERSION, 'label':label,
            'source_detection':{'predictor':s.get('predictor'),'version':s.get('predictor_version'),'model':s.get('predictor_model')},
            'translation_table_fallback':s['translation_table'],
            'translation_table_source':'explicit' if settings.get('translation_table') is not None else 'default_nuclear_assumption',
            'automatic_annevo_comment_lengths':s['automatic_annevo_lengths'],
            'output_directory':s['output_dir'],
            'policy':'Preserve predictions; add representation and submission metadata; biological QC is advisory.'}
        write_json(directory/'preparation.json',report)
        warnings = report['warnings']
        summary = [f'ENA annotation package {VERSION}',f'Label: {label}',
            f'Genes: {report["genes"]}',f'Feature rows: {report["feature_rows"]}',
            f'Inference: {s.get("inference") or "preserve existing only"}',
            f'Translation-table fallback: {s["translation_table"]}',
            f'Warnings: {len(warnings)}',
            'No genes or coordinates filtered for biological QC. Official Webin validation is still required.']
        summary += ['WARNING: '+w for w in warnings]
        (directory/'preparation.txt').write_text('\n'.join(summary)+'\n',encoding='ascii')
        validation_instructions(directory,'annotation')
    return {'status':'prepared','label':label,'output_directory':s['output_dir'],
            'genes':report['genes'],'feature_rows':report['feature_rows'],'warning_count':len(report['warnings'])}


def batch_prepare(path, defaults, summary_path):
    """Prepare TSV jobs sequentially and write a JSON outcome summary.
    
    Nonblank row values override defaults; row paths are relative to the TSV.
    Continue after individual preparation failures. Return 1 if any job failed,
    otherwise 0. Repeated prefixes are refused to avoid independent tag allocation."""
    rows = tsv_rows(path, SETTING_KEYS-{'output','report','manifest','write_locus_tag_map'}, ['input'])
    check(rows,'Batch TSV has no rows')
    summary = Path(summary_path)
    check(not summary.exists(),f'Batch summary already exists: {summary}')
    summary.parent.mkdir(parents=True,exist_ok=True)
    outcomes=[]; seen_dirs=set(); seen_prefixes=set()
    for number, row in enumerate(rows,2):
        settings = dict(defaults)
        settings.update(relative_paths({k:v for k,v in row.items() if v},Path(path).resolve().parent))
        # TSV values are strings; boolean fields have an explicit small vocabulary.
        try:
            for key in ['no_inference','strict_warnings','use_annevo_lengths']:
                if isinstance(settings.get(key),str):
                    check(settings[key].lower() in {'true','false'},f'{key} must be true or false')
                    settings[key]=settings[key].lower()=='true'
            resolved=resolve_settings(settings)
            check(resolved['output_dir'] not in seen_dirs,'Repeated batch output directory')
            check(resolved['locus_tag_prefix'] not in seen_prefixes,
                  'Repeated locus-tag prefix in this batch could allocate duplicate tags; use distinct registered prefixes or prepare coordinated updates separately')
            seen_dirs.add(resolved['output_dir']);seen_prefixes.add(resolved['locus_tag_prefix'])
            result=prepare_package(settings)
        except (ValueError,OSError,UnicodeError) as error:
            result={'status':'failed','error':str(error)}
        outcomes.append({'row':number,'input':settings.get('input'),**result})
        print(json.dumps(outcomes[-1]),flush=True)
    write_json(summary,{'toolkit_version':VERSION,'rows':outcomes,'prepared':sum(r['status']=='prepared' for r in outcomes)})
    return 1 if any(r['status']=='failed' for r in outcomes) else 0


def parser():
    """Build the CLI parser for package, batch, audit and legacy output modes."""
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--version',action='version',version=VERSION)
    for key in sorted(SETTING_KEYS-{'no_inference','strict_warnings','use_annevo_lengths','taxon_id','translation_table','width','start','step'}):
        p.add_argument('--'+key.replace('_','-'),
            help='Model used by the predictor; included in default label, manifest name and description.' if key=='predictor_model' else None)
    for key in ['taxon_id','translation_table','width','start','step']:
        p.add_argument('--'+key.replace('_','-'),type=int,
            help='Optional taxonomy ID; adds ##species. Without it, an existing valid directive is preserved.' if key=='taxon_id' else None)
    for key in ['no_inference','strict_warnings','use_annevo_lengths']:
        p.add_argument('--'+key.replace('_','-'),action='store_true',default=None)
    p.add_argument('--config',help='JSON defaults; paths resolve relative to the JSON file')
    p.add_argument('--batch',help='TSV of annotations; row paths resolve relative to the TSV')
    p.add_argument('--batch-summary',default='batch_summary.json')
    p.add_argument('--audit-only',action='store_true',help='Descriptive audit only; requires --input and --report')
    return p


def main(argv=None):
    """Dispatch CLI arguments and return an exit status.
    
    Return 0 on success, 1 for batch failures and 2 for caught preparation errors.
    CLI values override JSON defaults; batch rows can override those merged values."""
    p=parser();args=p.parse_args(argv)
    try:
        settings=configuration(args.config)
        settings.update({k:v for k,v in vars(args).items() if k in SETTING_KEYS and v is not None})
        if args.audit_only:
            check(not args.batch,'Audit mode does not accept --batch')
            check(settings.get('input') and settings.get('report'),'Audit requires --input and --report')
            check(Path(settings['input']).resolve()!=Path(settings['report']).resolve(),'Report cannot overwrite input')
            with core.staged_outputs([settings['report']]) as staged:
                write_json(staged[settings['report']],core.audit(settings['input']))
            return 0
        if args.batch:
            check(not any(settings.get(k) for k in ['output','report','manifest','write_locus_tag_map']),
                  'Legacy output-path options cannot be combined with --batch')
            return batch_prepare(args.batch,settings,args.batch_summary)
        legacy=['output','report','manifest','write_locus_tag_map']
        if any(settings.get(k) for k in legacy):
            check(all(settings.get(k) for k in ['output','report','write_locus_tag_map']),
                  'Legacy mode requires --output, --report and --write-locus-tag-map together; otherwise use --output-dir')
            s=resolve_settings(settings)
            for key in ['fai','assembly_report','locus_tag_map','manifest']:
                s.setdefault(key,None)
            report=core.prepare(SimpleNamespace(**s))
            print(json.dumps({'status':'prepared','genes':report['genes'],'output':s['output']}))
        else:
            print(json.dumps(prepare_package(settings),indent=2))
        return 0
    except (ValueError,OSError,UnicodeError,TypeError) as error:
        print('ERROR: '+str(error),file=sys.stderr);return 2

if __name__=='__main__':
    sys.exit(main())
