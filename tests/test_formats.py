"""Cross-dialect regression checks for graph, provenance and coding preservation."""
import gzip
import json
import tempfile
import unittest
from pathlib import Path
import ena_gffprep as g
from ena_formats import gtf_attributes

class FormatTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.input=self.root/'input.gff'
        self.settings=dict(input=str(self.input),assembly_accession='GCA_123456789.1',bioproject='PRJEB123456',locus_tag_prefix='TEST',no_inference=True,output_dir=str(self.root/'out'))
    def run_text(self, text, **kw):
        self.input.write_text('##gff-version 3\n##sequence-region s 1 1000\n'+text)
        result=g.prepare_package({**self.settings,**kw})
        p=Path(result['output_directory'])
        return gzip.open(next(p.glob('*.gff3.gz')),'rt').read(),json.loads((p/'preparation.json').read_text())
    def fixture(self, source='Tiberius'):
        return '\n'.join('\t'.join(x) for x in [
            ['s',source,'gene','1','100','.','+','.','ID=g'],
            ['s',source,'transcript','1','100','.','+','.','ID=t;Parent=g'],
            ['s',source,'CDS','1','9','.','+','0','ID=a;Parent=t'],
            ['s',source,'CDS','20','28','.','+','0','ID=b;Parent=t']])+'\n'
    def test_tiberius_coordinates_phases_sources_preserved(self):
        out,r=self.run_text(self.fixture())
        self.assertIn('\tmRNA\t',out);self.assertEqual(out.count('ID=t.CDS;'),2)
        self.assertIn('original_id=a',out);self.assertIn('original_id=b',out)
        self.assertIn('Tiberius\tCDS\t20\t28\t.\t+\t0',out)
    def test_unknown_version_not_invented(self):
        self.input.write_text(self.fixture())
        settings=g.resolve_settings({**self.settings,'no_inference':False})
        self.assertNotIn('predictor_version',settings)
        self.assertEqual(settings['inference_by_source']['Tiberius'],'ab initio prediction:Tiberius')
    def test_mixed_requires_inference_choice(self):
        self.input.write_text(self.fixture('Tiberius,transdecoder'))
        # Combined source tokens are handled as mixed evidence signatures.
        with self.assertRaises(ValueError):g.resolve_settings({**self.settings,'no_inference':False,'predictor':'Tiberius','predictor_version':'1','prediction_method':'evidence-assisted'})
    def test_braker_phases_and_transcript_attributes(self):
        text=self.fixture('AUGUSTUS')+'s\tAGAT\texon\t1\t9\t.\t+\t2\tID=e;Parent=t;transcript_id=t\n'
        out,r=self.run_text(text,predictor='BRAKER',predictor_version='4')
        self.assertIn('\texon\t1\t9\t.\t+\t.\t',out)
        self.assertIn('original_phase=2',out);self.assertIn('original_transcript_id=t',out)
        self.assertNotIn('Ab initio', (self.root/'out/manifest.txt').read_text())
    def test_gtf_quoted_metadata_and_graph(self):
        text='s\tOrionGeno\tgene\t1\t100\t.\t+\t.\tgene_id "g";\n'
        text+='s\tOrionGeno\ttranscript\t1\t100\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
        text+='s\tOrionGeno\tCDS\t1\t9\t.\t+\t0\tgene_id "g"; transcript_id "t"; note "a;b";\n'
        out,r=self.run_text(text)
        self.assertIn('Parent=t',out);self.assertIn('note=a%3Bb',out)
        with self.assertRaises(ValueError):gtf_attributes('gene_id g;')
    def test_conflicting_products_refused(self):
        with self.assertRaisesRegex(ValueError,'Conflicting product'):
            self.run_text(self.fixture().replace('ID=a;Parent=t','ID=a;Parent=t;product=A').replace('ID=b;Parent=t','ID=b;Parent=t;product=B'))
    def test_carbon_grouped_cds_retained(self):
        out,r=self.run_text(self.fixture('Carbon-A').replace('ID=a;','ID=c;').replace('ID=b;','ID=c;'))
        self.assertEqual(out.count('ID=c;'),2)
    def test_helixer_name_not_guessed(self):
        out,r=self.run_text('##species Conistra_vaccinii\n'+self.fixture('Helixer'))
        self.assertNotIn('##species',out)
        self.assertTrue(any('Name-only' in w for w in r['warnings']))
    def test_isoforms_remain_separate(self):
        text=self.fixture()+self.fixture().split('\n',1)[1].replace('ID=t;','ID=u;').replace('Parent=t','Parent=u').replace('ID=a;','ID=d;').replace('ID=b;','ID=e;')
        out,r=self.run_text(text)
        self.assertEqual(r['genes'],1);self.assertIn('ID=t.CDS;',out);self.assertIn('ID=u.CDS;',out)
    def test_source_map_applies_only_to_selected_sources(self):
        mapping=self.root/'inference.json';mapping.write_text(json.dumps({'Tiberius':'ab initio prediction:Tiberius:1.0'}))
        out,r=self.run_text(self.fixture(),no_inference=False,inference_map=str(mapping))
        self.assertEqual(out.count('inference=ab initio prediction:Tiberius:1.0'),2)
    def test_source_map_unknown_source_refused(self):
        mapping=self.root/'inference.json';mapping.write_text('{"typo":"ab initio prediction:X:1"}')
        with self.assertRaisesRegex(ValueError,'unobserved'):self.run_text(self.fixture(),no_inference=False,inference_map=str(mapping))
    def test_gzip_reproducible_across_paths(self):
        self.run_text(self.fixture())
        a=next((self.root/'out').glob('*.gz')).read_bytes()
        g.prepare_package({**self.settings,'output_dir':str(self.root/'second')})
        b=next((self.root/'second').glob('*.gz')).read_bytes()
        self.assertEqual(a,b)
    def test_conflicting_redundant_identifier_refused(self):
        text=self.fixture()+'s\tTiberius\texon\t1\t9\t.\t+\t.\tID=e;Parent=t;transcript_id=other\n'
        with self.assertRaisesRegex(ValueError,'transcript_id conflicts'):self.run_text(text)
    def test_generic_remains_strict(self):
        text=self.fixture()+'s\tAUGUSTUS\texon\t1\t9\t.\t+\t1\tID=e;Parent=t\n'
        with self.assertRaisesRegex(ValueError,'invalid phase'):self.run_text(text,source_format='generic')
    def test_combined_source_detected_as_assisted(self):
        self.input.write_text(self.fixture('Tiberius,transdecoder'))
        s=g.resolve_settings(self.settings)
        self.assertEqual(s['source_format'],'tiberius');self.assertEqual(s['prediction_method'],'evidence-assisted')
    def test_explicit_taxonomy_replaces_helixer_name(self):
        out,r=self.run_text('##species Conistra_vaccinii\n'+self.fixture('Helixer'),taxon_id=123)
        self.assertIn('##species https://www.ncbi.nlm.nih.gov/datasets/taxonomy/123/',out)
    def test_split_codons_and_minus_strand_preserved(self):
        text=self.fixture('AUGUSTUS').replace('\t+\t','\t-\t')
        text+='s\tAUGUSTUS\tstart_codon\t27\t28\t.\t-\t0\tID=start1;Parent=t\n'
        text+='s\tAUGUSTUS\tstart_codon\t9\t9\t.\t-\t0\tID=start2;Parent=t\n'
        out,r=self.run_text(text)
        self.assertIn('\tstart_codon\t27\t28\t.\t-\t.\t',out)
        self.assertIn('\tCDS\t20\t28\t.\t-\t0',out)
    def test_idempotent_representation(self):
        out,r=self.run_text(self.fixture())
        self.input.write_text(out)
        g.prepare_package({**self.settings,'output_dir':str(self.root/'second')})
        self.assertEqual(out,gzip.open(next((self.root/'second').glob('*.gz')),'rt').read())
    def test_stop_codon_non_cds_phase_normalized(self):
        out,r=self.run_text(self.fixture('AUGUSTUS')+'s\tAUGUSTUS\tstop_codon\t26\t28\t.\t+\t0\tID=stop;Parent=t\n')
        self.assertIn('\tstop_codon\t26\t28\t.\t+\t.\t',out)
    def test_map_relative_config_path(self):
        mapping=self.root/'evidence.json';mapping.write_text('{"Tiberius":null}')
        config=self.root/'defaults.json';config.write_text('{"inference_map":"evidence.json"}')
        self.assertEqual(g.configuration(config)['inference_map'],str(mapping))
    def test_map_conflicting_existing_evidence_refused(self):
        mapping=self.root/'evidence.json';mapping.write_text('{"Tiberius":"ab initio prediction:Tiberius:1.0"}')
        with self.assertRaisesRegex(ValueError,'existing inference differs'):
            self.run_text(self.fixture().replace('Parent=t','Parent=t;inference=ab initio prediction:Other:1'),no_inference=False,inference_map=str(mapping))
    def test_mixed_auto_rejects_blanket_inference(self):
        self.input.write_text(self.fixture('Tiberius,transdecoder'))
        with self.assertRaisesRegex(ValueError,'explicit --inference'):
            g.resolve_settings({**self.settings,'no_inference':False,'predictor':'Tiberius','predictor_version':'1.0'})
    def test_direct_sources_in_mixed_input_get_ab_initio_only(self):
        text=self.fixture()+self.fixture('transdecoder').split('\n',1)[1].replace('ID=t;','ID=u;').replace('Parent=t','Parent=u').replace('ID=a;','ID=d;').replace('ID=b;','ID=e;')
        out,r=self.run_text(text,no_inference=False,predictor_version='2.0')
        rows=[line.split('\t') for line in out.splitlines() if line and not line.startswith('#')]
        direct=[x for x in rows if x[1]=='Tiberius' and x[2]=='CDS']
        other=[x for x in rows if x[1]=='transdecoder' and x[2]=='CDS']
        self.assertTrue(all('inference=ab initio prediction:Tiberius:2.0' in x[8] for x in direct))
        self.assertTrue(all('inference=' not in x[8] for x in other))
    def test_pipeline_version_not_reused_as_component_version(self):
        self.input.write_text(self.fixture())
        settings=g.resolve_settings({**self.settings,'no_inference':False,'predictor':'MyPipeline','predictor_version':'9'})
        self.assertEqual(settings['inference_by_source']['Tiberius'],'ab initio prediction:Tiberius')
    def test_automatic_source_policy_preserves_existing_inference(self):
        out,r=self.run_text(self.fixture().replace('Parent=t','Parent=t;inference=similar to RNA sequence:INSD:ABC123'),no_inference=False)
        self.assertIn('inference=similar to RNA sequence:INSD:ABC123',out)
        self.assertNotIn('inference=ab initio',out)
