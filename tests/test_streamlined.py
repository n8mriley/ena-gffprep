import contextlib
import gzip
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import ena_gffprep as g

GENE='''##gff-version 3
##sequence-region seq1 1 100
seq1\ttest\tgene\t1\t90\t.\t+\t.\tID=g1
seq1\ttest\tmRNA\t1\t90\t.\t+\t.\tID=t1;Parent=g1
seq1\ttest\tCDS\t1\t90\t.\t+\t0\tID=c1;Parent=t1
'''
ANNEVO=GENE.replace('ID=c1','ID=t1.CDS.1')
ANNEVO='# This output was generated with ANNEVO v2.2.3.\n'+ANNEVO

class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        (self.root/'source.gff').write_text(ANNEVO)
        self.settings={'input':str(self.root/'source.gff'),'assembly_accession':'GCA_123456789.1',
            'taxon_id':9606,'locus_tag_prefix':'TEST','output_dir':str(self.root/'prepared')}
    def test_auto_package_and_manifest(self):
        # This fixture uses directives, not ANNEVO length comments: supply an external index.
        (self.root/'a.fai').write_text('seq1\t100\n')
        r=g.prepare_package({**self.settings,'fai':str(self.root/'a.fai')})
        p=self.root/'prepared';self.assertEqual(r['genes'],1)
        self.assertTrue((p/'manifest.txt').is_file());self.assertTrue((p/'VALIDATE.txt').is_file())
        report=json.loads((p/'preparation.json').read_text())
        self.assertEqual(report['workflow']['source_detection']['version'],'2.2.3')
        target=[x for x in (p/'manifest.txt').read_text().splitlines() if x.startswith('GFF3\t')][0].split('\t')[1]
        with gzip.open(p/target,'rt') as h:out=h.read()
        self.assertIn('ID=t1.CDS;',out);self.assertIn('inference=ab initio prediction:ANNEVO:2.2.3',out)
        with self.assertRaisesRegex(ValueError,'already exists'):g.prepare_package({**self.settings,'fai':str(self.root/'a.fai')})
    def test_provenance_conflict(self):
        with self.assertRaisesRegex(ValueError,'version conflicts'):g.resolve_settings({**self.settings,'predictor_version':'2.2.1'})
        with self.assertRaisesRegex(ValueError,'inference conflicts'):g.resolve_settings({**self.settings,'inference':'ab initio prediction:Other:1'})
    def test_generic_requires_explicit_provenance(self):
        (self.root/'generic.gff').write_text(GENE)
        s={**self.settings,'input':str(self.root/'generic.gff')}
        with self.assertRaisesRegex(ValueError,'Cannot determine'):g.resolve_settings(s)
        r=g.resolve_settings({**s,'predictor':'Predictor','predictor_version':'1.0'})
        self.assertEqual(r['source_format'],'generic')
        self.assertEqual(r['inference'],'ab initio prediction:Predictor:1.0')
    def test_failed_package_is_not_published(self):
        (self.root/'source.gff').write_text(ANNEVO.replace('Parent=t1','Parent=missing'))
        (self.root/'a.fai').write_text('seq1\t100\n')
        with self.assertRaises(ValueError):g.prepare_package({**self.settings,'fai':str(self.root/'a.fai')})
        self.assertFalse((self.root/'prepared').exists());self.assertFalse(list(self.root.glob('.ena-prep-*')))
    def test_batch_continues_and_resolves_paths(self):
        (self.root/'a.fai').write_text('seq1\t100\n')
        (self.root/'batch.tsv').write_text('input\tfai\tlocus_tag_prefix\tlabel\nsource.gff\ta.fai\tTESA\ta\nmissing.gff\ta.fai\tTESB\tb\nsource.gff\ta.fai\tTESC\tc\n')
        defaults={k:v for k,v in self.settings.items() if k not in {'input','output_dir','locus_tag_prefix'}}
        defaults['out_root']=str(self.root/'batch')
        with contextlib.redirect_stdout(io.StringIO()):code=g.batch_prepare(self.root/'batch.tsv',defaults,self.root/'batch.json')
        self.assertEqual(code,1)
        results=json.loads((self.root/'batch.json').read_text());self.assertEqual(results['prepared'],2)
        self.assertEqual([r['status'] for r in results['rows']],['prepared','failed','prepared'])
    def test_config_paths_and_unknown_settings(self):
        p=self.root/'config.json';p.write_text(json.dumps({'fai':'a.fai','translation_table':1}))
        self.assertEqual(g.configuration(p)['fai'],str(self.root/'a.fai'))
        p.write_text('{"typo":1}')
        with self.assertRaises(ValueError):g.configuration(p)


class AdditionalSafeguards(unittest.TestCase):
    def test_annevo_comment_lengths_without_index(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);source=p/'source.gff'
            text=ANNEVO.replace('##sequence-region seq1 1 100\n','# ----- prediction on sequence number 1 (length = 100, name = seq1) -----\n')
            source.write_text(text)
            r=g.prepare_package({'input':str(source),'assembly_accession':'GCA_123456789.1',
                'taxon_id':9606,'locus_tag_prefix':'TEST','output_dir':str(p/'out')})
            self.assertEqual(r['genes'],1)
            report=json.loads((p/'out/preparation.json').read_text())
            self.assertTrue(report['workflow']['automatic_annevo_comment_lengths'])
            self.assertTrue(any('only' in w and 'length' in w.lower() for w in report['warnings']))
