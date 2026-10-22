"""Offline organelle filtering, network failure handling and publication checks."""
import contextlib
import copy
import gzip
import io
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from urllib.error import HTTPError, URLError
import ena_organelles as o
import ena_gffprep as cli
import test_gffprep as base

ACC='GCA_123456789.1'

def record(gb, location='Chromosome', name=None, rs=None, **extra):
    """Create a synthetic NCBI sequence record, never querying the network."""
    return dict(assembly_accession=ACC, genbank_accession=gb, sequence_name=name or gb,
                refseq_accession=rs, assigned_molecule_location_type=location,
                assembly_unit='Primary Assembly' if location=='Chromosome' else 'non-nuclear',
                length=1000, role='assembled-molecule', **extra)

MT='mt\ttest\tgene\t1\t90\t.\t+\t.\tID=gm\nmt\ttest\tmRNA\t1\t90\t.\t+\t.\tID=tm;Parent=gm\nmt\ttest\tCDS\t1\t90\t.\t+\t0\tID=cm;Parent=tm\n'
TEXT=base.BASE.replace('##sequence-region seq1 1 1000','##sequence-region seq1 1 1000\n##sequence-region mt 1 1000')+MT

class OrganelleTests(unittest.TestCase):
    setUp=base.Tests.setUp
    def runprep(self, text=base.BASE, **kwargs):
        kwargs.setdefault("assembly_accession",ACC)
        return base.Tests.runprep(self,text,**kwargs)

    def metadata(self, records=None):
        p=self.d/'sequence.json'
        p.write_text(json.dumps({'reports':records if records is not None else [record('seq1'),record('mt','Mitochondrion')]}))
        return str(p)

    def prep(self,text=TEXT,**extra):
        return self.runprep(text,remove_organelles=True,organelle_metadata=self.metadata(),**extra)

    def test_mixed_filter_and_headers(self):
        a,r,out=self.prep()
        self.assertNotIn('##sequence-region mt ',out)
        self.assertNotIn('\nmt\t',out)
        self.assertIn('##sequence-region seq1 ',out)
        self.assertEqual(r['organelle_filter']['removed_feature_rows'],3)
        self.assertEqual(r['organelle_filter']['removed_genes'],1)
        report=Path(a.report).with_suffix('.organelles.tsv').read_text()
        self.assertIn('mt\tmt\tMitochondrion\t3\t1',report)

    def test_flag_off_no_lookup_and_identical_gff(self):
        _,_,before=self.runprep(TEXT)
        with patch.object(o,'load_metadata',side_effect=AssertionError('Network must not be used')):
            _,r,after=self.runprep(TEXT,remove_organelles=False)
        self.assertEqual(before,after)
        self.assertNotIn('organelle_filter',r)

    def test_alias_matching_refseq_and_original(self):
        records=[record('GB1.1',name='seq1',rs='RS1.1'),record('GBMT.1','Mitochondrion',name='mt',rs='RSMT.1')]
        _,r,out=self.runprep(TEXT.replace('mt','RSMT.1'),remove_organelles=True,organelle_metadata=self.metadata(records))
        self.assertNotIn('RSMT.1',out)
        self.assertEqual(r['organelle_filter']['excluded_sequences'][0]['ncbi_sequence'],'GBMT.1')

    def test_no_organelle_results(self):
        a,r,out=self.runprep(base.BASE,remove_organelles=True,organelle_metadata=self.metadata([record('seq1')]))
        self.assertEqual(r['organelle_filter']['removed_feature_rows'],0)
        self.assertEqual(len(Path(a.report).with_suffix('.organelles.tsv').read_text().splitlines()),1)

    def test_organelle_in_assembly_but_unannotated(self):
        _,r,_=self.prep(base.BASE)
        self.assertEqual(r['organelle_filter']['removed_feature_rows'],0)
        self.assertEqual(len(r['organelle_filter']['organelles_in_assembly']),1)

    def test_plastids_and_other_compartments(self):
        for location in ['Chloroplast','Plastid','Apicoplast','Kinetoplast','Plastid:chloroplast']:
            with self.subTest(location=location):
                _,r,_=self.runprep(TEXT,remove_organelles=True,organelle_metadata=self.metadata([record('seq1'),record('mt',location)]))
                self.assertEqual(r['organelle_filter']['removed_feature_rows'],3)

    def test_wrong_assembly_and_version(self):
        for wrong in ['GCA_123456789.2','GCF_123456789.1']:
            records=[record('seq1')];records[0]['assembly_accession']=wrong
            with self.assertRaisesRegex(ValueError,'accession/version'):
                self.runprep(base.BASE,remove_organelles=True,organelle_metadata=self.metadata(records))

    def test_unmatched_and_version_not_stripped(self):
        for name in ['missing','seq1.2']:
            with self.assertRaisesRegex(ValueError,'Unmatched'):
                self.runprep(base.BASE.replace('seq1',name),remove_organelles=True,organelle_metadata=self.metadata())

    def test_ambiguous_alias(self):
        records=[record('seq1'),record('other','Mitochondrion',name='seq1')]
        with self.assertRaisesRegex(ValueError,'Ambiguous'):
            self.runprep(base.BASE,remove_organelles=True,organelle_metadata=self.metadata(records))

    def test_all_removed(self):
        with self.assertRaisesRegex(ValueError,'remove all'):
            self.runprep('##gff-version 3\n##sequence-region mt 1 1000\n'+MT,remove_organelles=True,organelle_metadata=self.metadata())

    def test_cross_boundary_parent_and_derives_from(self):
        for attr in ['Parent=gm','Derives_from=gm']:
            with self.subTest(attr=attr), self.assertRaisesRegex(ValueError,'references an excluded'):
                self.prep(TEXT.replace('ID=g1','ID=g1;'+attr))

    def test_shared_id_across_boundary(self):
        with self.assertRaisesRegex(ValueError,'shared by retained'):
            self.prep(TEXT.replace('ID=gm','ID=g1'))

    def test_existing_excluded_tags_reserved_without_map(self):
        a,r,out=self.prep(TEXT.replace('ID=gm','ID=gm;locus_tag=TEST_000010'))
        self.assertIn('locus_tag=TEST_000011',out)
        self.assertIn('gm\tTEST_000010\treserved',Path(a.write_locus_tag_map).read_text())

    def test_excluded_child_tag_reserved(self):
        a,_,_=self.prep(TEXT.replace('ID=cm;Parent=tm','ID=cm;Parent=tm;locus_tag=TEST_A1'))
        self.assertIn('gm\tTEST_A1\treserved',Path(a.write_locus_tag_map).read_text())

    def test_previous_map_retains_removed_and_active_tags(self):
        p=self.d/'previous.tsv';p.write_text('gene_id\tlocus_tag\ng1\tTEST_000003\ngm\tTEST_000009\nold\tTEST_Z1\n')
        a,_,out=self.prep(locus_tag_map=str(p))
        self.assertIn('locus_tag=TEST_000003',out)
        text=Path(a.write_locus_tag_map).read_text()
        self.assertIn('gm\tTEST_000009\treserved',text)
        self.assertIn('old\tTEST_Z1\treserved',text)

    def test_missing_or_unknown_classification_fails(self):
        for location in ['', 'Imaginary organelle']:
            with self.assertRaises(ValueError):o.validate_records([record('seq1',location)],ACC)

    def test_unclassified_primary_scaffold_reported(self):
        r=record('seq1');r.pop('assigned_molecule_location_type');r['role']='unplaced-scaffold'
        _,report,out=self.runprep(base.BASE,remove_organelles=True,organelle_metadata=self.metadata([r]))
        self.assertEqual(report['organelle_filter']['unclassified_retained_sequences'],['seq1'])
        self.assertTrue(any('unspecified' in w for w in report['warnings']))

    def test_length_mismatch_fails(self):
        records=[record('seq1')];records[0]['length']=2000
        with self.assertRaisesRegex(ValueError,'length mismatch'):
            self.runprep(base.BASE,remove_organelles=True,organelle_metadata=self.metadata(records))

    def test_jsonl_camel_case(self):
        p=self.d/'sequences.jsonl';p.write_text(json.dumps(dict(assemblyAccession=ACC,genbankAccession='seq1',length=1000,assignedMoleculeLocationType='Chromosome'))+'\n')
        _,r,_=self.runprep(base.BASE,remove_organelles=True,organelle_metadata=str(p))
        self.assertEqual(r['organelle_filter']['metadata_sequence_count'],1)

    def test_partial_and_empty_offline_report_fails(self):
        for data in [{'reports':[]},{'reports':[record('seq1')],'next_page_token':'more'},{'reports':[record('seq1')],'total_count':2}]:
            p=self.d/'bad.json';p.write_text(json.dumps(data))
            with self.assertRaises(ValueError):self.runprep(base.BASE,remove_organelles=True,organelle_metadata=str(p))

    def test_package_failure_not_published(self):
        dest=self.d/'package'
        source=self.d/'input.gff';source.write_text(base.BASE)
        with self.assertRaisesRegex(ValueError,'Unmatched'):
            cli.prepare_package(dict(input=str(source),organelle_metadata=self.metadata([record('mt','Mitochondrion')]),assembly_accession=ACC,locus_tag_prefix='TEST',no_inference=True,remove_organelles=True,output_dir=str(dest)))
        self.assertFalse(dest.exists())

class MetadataTransportTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.args=SimpleNamespace(assembly_accession=ACC,organelle_cache_dir=str(self.root/'cache'))

    def test_fetch_pagination_and_cache(self):
        pages=[{'reports':[record('seq1')],'total_count':2,'next_page_token':'next'},
               {'reports':[record('mt','Mitochondrion')],'total_count':2}]
        with patch.object(o,'request_json',side_effect=pages) as request:
            records,meta=o.load_metadata(self.args)
            self.assertEqual(request.call_count,2)
            self.assertIn('page_token=next',request.call_args[0][0])
        self.assertEqual(meta['source'],'ncbi_api');self.assertEqual(len(records),2)
        self.assertTrue(meta['retrieved_at']);self.assertEqual(len(meta['metadata_sha256']),64)
        with patch.object(o,'request_json',side_effect=AssertionError('must use cache')):
            _,again=o.load_metadata(self.args)
        self.assertEqual(again['source'],'cache')

    def test_refresh_cache(self):
        with patch.object(o,'request_json',return_value={'reports':[record('seq1')]}):o.load_metadata(self.args)
        self.args.refresh_organelle_cache=True
        with patch.object(o,'request_json',return_value={'reports':[record('mt','Mitochondrion')]}) as request:
            records,_=o.load_metadata(self.args)
        self.assertTrue(request.called);self.assertTrue(records[0]['organelle'])

    def test_incomplete_or_repeating_network_pages(self):
        cases=[[{'reports':[record('seq1')],'total_count':2}],
               [{'reports':[record('seq1')],'next_page_token':'x'}, {'reports':[],'next_page_token':'x'}]]
        for pages in cases:
            with self.subTest(pages=pages), patch.object(o,'request_json',side_effect=pages), self.assertRaises(ValueError):
                o.load_metadata(self.args)
        self.assertFalse((self.root/'cache'/f'{ACC}.json').exists())

    def test_retry_transient_failure(self):
        with patch.object(o,'urlopen',side_effect=[URLError('offline'),io.StringIO(json.dumps({'reports':[]}))]) as request, patch.object(o.time,'sleep'):
            o.request_json('https://example.invalid')
        self.assertEqual(request.call_count,2)

    def test_network_failure_no_cache(self):
        with patch.object(o,'urlopen',side_effect=URLError('offline')), patch.object(o.time,'sleep'), self.assertRaisesRegex(ValueError,'NCBI metadata request failed'):
            o.load_metadata(self.args)
        self.assertFalse((self.root/'cache'/f'{ACC}.json').exists())

    def test_permanent_http_failure_not_retried(self):
        error=HTTPError('https://example.invalid',404,'not found',{},None)
        with patch.object(o,'urlopen',side_effect=error) as request, self.assertRaisesRegex(ValueError,'HTTP 404'):
            o.load_metadata(self.args)
        self.assertEqual(request.call_count,1)

    def test_corrupt_cache_fails_without_network(self):
        path=self.root/'cache';path.mkdir();(path/f'{ACC}.json').write_text('not JSON')
        with patch.object(o,'request_json',side_effect=AssertionError('must not silently replace')), self.assertRaises(ValueError):
            o.load_metadata(self.args)

    def test_cli_modes_and_config_paths(self):
        a=cli.parser().parse_args(['--remove-organelles','--organelle-metadata','sequence.json'])
        self.assertTrue(a.remove_organelles);self.assertEqual(a.organelle_metadata,'sequence.json')
        p=self.root/'settings.json';p.write_text(json.dumps({'organelle_metadata':'sequence.json','remove_organelles':True}))
        self.assertEqual(cli.configuration(p)['organelle_metadata'],str(self.root/'sequence.json'))
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(cli.main(['--audit-only','--remove-organelles']),2)

    def test_batch_filter_and_sidecar(self):
        metadata=self.root/'sequence.json';metadata.write_text(json.dumps({'reports':[record('seq1'),record('mt','Mitochondrion')]}))
        source=self.root/'input.gff';source.write_text(TEXT)
        batch=self.root/'jobs.tsv'
        batch.write_text('input\tassembly_accession\tlocus_tag_prefix\tno_inference\tremove_organelles\torganelle_metadata\ninput.gff\t'+ACC+'\tTEST\ttrue\ttrue\tsequence.json\n')
        with contextlib.redirect_stdout(io.StringIO()):
            code=cli.batch_prepare(batch,{'out_root':str(self.root/'prepared')},self.root/'summary.json')
        self.assertEqual(code,0)
        outputs=list((self.root/'prepared').glob('*/organelle_exclusions.tsv'))
        self.assertEqual(len(outputs),1)
        self.assertIn('Mitochondrion',outputs[0].read_text())

    def test_options_without_flag_rejected(self):
        source=self.root/'input.gff';source.write_text(base.BASE)
        with self.assertRaisesRegex(ValueError,'require --remove-organelles'):
            cli.resolve_settings(dict(input=str(source),assembly_accession=ACC,locus_tag_prefix='TEST',no_inference=True,organelle_metadata='meta.json'))
