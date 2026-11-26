"""Length fallback, source precedence and shared metadata without network access."""
import contextlib
import gzip
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import ena_gffprep as g
import ena_organelles as ena
from test_gffprep import BASE
from test_organelles import ACC, record, TEXT, tsv

class LengthFallbackTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.input=self.root/'input.gff'
        self.settings=dict(input=str(self.input),assembly_accession=ACC,bioproject='PRJEB123456',locus_tag_prefix='TEST',no_inference=True,output_dir=str(self.root/'out'),sequence_cache_dir=str(self.root/'cache'))
        self.missing='\n'.join(l for l in BASE.splitlines() if not l.startswith('##sequence-region'))+'\n'
    def prep(self,text=None,**extra):
        self.input.write_text(self.missing if text is None else text)
        result=g.prepare_package({**self.settings,**extra})
        dest=Path(result['output_directory'])
        return json.loads((dest/'preparation.json').read_text()),gzip.open(next(dest.glob('*.gz')),'rt').read()
    def offline(self,records):
        p=self.root/(ACC+'_sequence_report.txt');p.write_bytes(tsv(records));return str(p)
    def test_automatic_lookup_and_provenance(self):
        with patch.object(ena,'request_report',return_value=tsv([record('seq1')])) as request:
            report,out=self.prep()
        self.assertEqual(request.call_count,1)
        self.assertEqual(report['sequence_lengths_source'],'ena_sequence_report')
        self.assertEqual(len(report['ena_sequence_lengths']['metadata_sha256']),64)
        self.assertIn('##sequence-region seq1.1 1 1000',out)
    def test_complete_directives_skip_lookup(self):
        with patch.object(ena,'load_metadata',side_effect=AssertionError('No lookup needed')):
            report,_=self.prep(BASE)
        self.assertNotIn('ena_sequence_lengths',report)
    def test_explicit_fai_does_not_trigger_lookup(self):
        p=self.root/'a.fai';p.write_text('seq1\t1000\n')
        with patch.object(ena,'load_metadata',side_effect=AssertionError('No lookup needed')):
            report,_=self.prep(fai=str(p))
        self.assertEqual(report['sequence_lengths_source'],'fai')
    def test_incomplete_fai_is_not_hidden(self):
        p=self.root/'a.fai';p.write_text('other\t1000\n')
        with patch.object(ena,'load_metadata',side_effect=AssertionError('No lookup allowed')):
            with self.assertRaisesRegex(ValueError,'missing from FAI'):self.prep(fai=str(p))
    def test_partial_directives_checked_against_fallback(self):
        second=self.missing.replace('seq1','seq2').replace('g1','g2').replace('t1','t2').replace('c1','c2')
        report,out=self.prep(BASE+second,sequence_metadata=self.offline([record('seq1'),record('seq2')]))
        self.assertEqual(report['genes'],2);self.assertIn('##sequence-region seq2.1 1 1000',out)
    def test_existing_length_disagreement_fails(self):
        second=self.missing.replace('seq1','seq2').replace('g1','g2').replace('t1','t2').replace('c1','c2')
        with self.assertRaisesRegex(ValueError,'length mismatch'):
            self.prep(BASE.replace('1 1000','1 1100')+second,sequence_metadata=self.offline([record('seq1'),record('seq2')]))
    def test_exact_version_required(self):
        p=self.root/'GCA_123456789.2_sequence_report.txt';p.write_bytes(tsv([record('seq1')]))
        with self.assertRaisesRegex(ValueError,'accession/version'):self.prep(sequence_metadata=str(p))
    def test_unmatched_alias_rejected(self):
        with self.assertRaisesRegex(ValueError,'Unmatched'):self.prep(sequence_metadata=self.offline([record('seq1.2',name='other')]))
    def test_ambiguous_alias_rejected(self):
        with self.assertRaisesRegex(ValueError,'Ambiguous'):
            self.prep(sequence_metadata=self.offline([record('GB1.1',name='seq1'),record('GB2.1',name='seq1')]))
    def test_alias_mapping_to_genbank(self):
        report,out=self.prep(sequence_metadata=self.offline([record('GB1.1',name='seq1',rs='RS1.1')]))
        self.assertIn('##sequence-region GB1.1 1 1000',out)
        self.assertIn('\nGB1.1\ttest\tgene',out)
        self.assertEqual(report['ena_sequence_lengths']['matched_sequences']['seq1'],'GB1.1')
    def test_refseq_column_rejected(self):
        with self.assertRaisesRegex(ValueError,'Only INSDC'):
            self.prep(accession_column='refseq')
    def test_feature_bounds_checked(self):
        short=record('seq1');short['sequence-length']='100'
        with self.assertRaisesRegex(ValueError,'beyond sequence length'):self.prep(sequence_metadata=self.offline([short]))
    def test_no_lookup_switch(self):
        with patch.object(ena,'load_metadata',side_effect=AssertionError('Offline policy')):
            with self.assertRaisesRegex(ValueError,'No sequence lengths'):self.prep(no_ena_lengths=True)
    def test_classification_not_needed_for_lengths(self):
        r=record('seq1');r['replicon-type']=''
        report,_=self.prep(sequence_metadata=self.offline([r]))
        self.assertEqual(report['genes'],1)
    def test_filtering_still_requires_classification(self):
        r=record('seq1');r['replicon-type']=''
        with self.assertRaisesRegex(ValueError,'classification'):
            self.prep(sequence_metadata=self.offline([r]),remove_organelles=True)
    def test_filter_and_lengths_share_one_request(self):
        text='\n'.join(l for l in TEXT.splitlines() if not l.startswith('##sequence-region'))+'\n'
        with patch.object(ena,'request_report',return_value=tsv([record('seq1'),record('mt','Mitochondrion')])) as request:
            report,out=self.prep(text,remove_organelles=True)
        self.assertEqual(request.call_count,1)
        self.assertTrue(report['ena_sequence_lengths']['reused_organelle_metadata'])
        self.assertEqual(report['organelle_filter']['removed_genes'],1)
        self.assertNotIn('\nmt\t',out)
    def test_cache_reused_between_packages(self):
        with patch.object(ena,'request_report',return_value=tsv([record('seq1')])):self.prep()
        with patch.object(ena,'request_report',side_effect=AssertionError('Use cache')):
            report,_=self.prep(output_dir=str(self.root/'second'))
        self.assertEqual(report['ena_sequence_lengths']['source'],'cache')
    def test_network_failure_does_not_publish(self):
        with patch.object(ena,'request_report',side_effect=ValueError('ENA metadata request failed')):
            with self.assertRaisesRegex(ValueError,'request failed'):self.prep()
        self.assertFalse((self.root/'out').exists())
    def test_annevo_missing_comments_falls_back(self):
        text=self.missing.replace('ID=c1','ID=t1.CDS')
        report,_=self.prep(text,predictor='ANNEVO',predictor_version='2.2.1',sequence_metadata=self.offline([record('seq1')]))
        self.assertEqual(report['sequence_lengths_source'],'ena_sequence_report')
    def test_annevo_comments_keep_existing_behavior(self):
        text=self.missing.replace('ID=c1','ID=t1.CDS')+'# ----- prediction on sequence number 1 (length = 1000, name = seq1) -----\n'
        with patch.object(ena,'load_metadata',side_effect=AssertionError('Comments suffice')):
            report,_=self.prep(text,predictor='ANNEVO',predictor_version='2.2.1')
        self.assertEqual(report['sequence_lengths_source'],'annevo_comments')
    def test_config_paths_and_boolean_flags(self):
        config=self.root/'defaults.json';config.write_text('{"sequence_metadata":"seq.json","sequence_cache_dir":"cache","no_ena_lengths":true}')
        settings=g.configuration(config)
        self.assertEqual(settings['sequence_metadata'],str(self.root/'seq.json'))
        args=g.parser().parse_args(['--no-ena-lengths','--refresh-sequence-cache'])
        self.assertTrue(args.no_ena_lengths);self.assertTrue(args.refresh_sequence_cache)
    def test_offline_refresh_conflict(self):
        self.input.write_text(self.missing)
        with self.assertRaisesRegex(ValueError,'Offline sequence metadata'):
            g.resolve_settings({**self.settings,'sequence_metadata':'seq.json','refresh_sequence_cache':True})
    def test_distinct_shared_metadata_rejected(self):
        self.input.write_text(self.missing)
        with self.assertRaisesRegex(ValueError,'one shared'):
            g.resolve_settings({**self.settings,'remove_organelles':True,'sequence_metadata':'seq.json','organelle_metadata':'other.json'})
    def test_batch_uses_offline_fallback(self):
        self.input.write_text(self.missing);self.offline([record('seq1')])
        batch=self.root/'jobs.tsv';batch.write_text('input\tassembly_accession\tlocus_tag_prefix\tno_inference\tsequence_metadata\ninput.gff\t'+ACC+'\tTEST\ttrue\t'+ACC+'_sequence_report.txt\n')
        with contextlib.redirect_stdout(io.StringIO()):
            code=g.batch_prepare(batch,{'bioproject':'PRJEB123456','out_root':str(self.root/'batch')},self.root/'summary.json')
        self.assertEqual(code,0)
