"""Predictor-led adapter selection, model provenance and ANNEVO header variants."""
import unittest
import tempfile
import json
from pathlib import Path
import ena_gffprep as g
from test_streamlined import ANNEVO

class PredictorModelTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.source = self.root/'input.gff'
        self.source.write_text(ANNEVO.replace('##sequence-region seq1 1 100', '# ----- prediction on sequence number 1 (length = 100, name = seq1) -----'))
        self.settings = dict(input=str(self.source), assembly_accession='GCA_123456789.1', locus_tag_prefix='TEST', output_dir=str(self.root/'out'))

    def test_unversioned_header_with_explicit_predictor(self):
        self.source.write_text(self.source.read_text().replace('ANNEVO v2.2.3.', 'ANNEVO.'))
        s = g.resolve_settings({**self.settings, 'predictor':'ANNEVO', 'predictor_version':'2.2.1'})
        self.assertEqual(s['source_format'], 'annevo')
        self.assertTrue(s['use_annevo_lengths'])
        g.prepare_package({**self.settings, 'predictor':'ANNEVO', 'predictor_version':'2.2.1'})
        self.assertTrue((self.root/'out/GCA_123456789.1.ANNEVO.2.2.1.ena.gff3.gz').exists())

    def test_explicit_predictor_without_header_selects_adapter(self):
        self.source.write_text('\n'.join(l for l in self.source.read_text().splitlines() if not l.startswith('# This output'))+'\n')
        s = g.resolve_settings({**self.settings, 'predictor':'annevo', 'predictor_version':'2.2.1'})
        self.assertEqual(s['source_format'], 'annevo')

    def test_unversioned_header_still_requires_version(self):
        self.source.write_text(self.source.read_text().replace('ANNEVO v2.2.3.', 'ANNEVO.'))
        with self.assertRaisesRegex(ValueError, 'Cannot determine'):
            g.resolve_settings(self.settings)

    def test_model_package_metadata(self):
        g.prepare_package({**self.settings, 'predictor_model':'Insecta'})
        label = 'GCA_123456789.1.ANNEVO.2.2.3.Insecta'
        self.assertTrue((self.root/'out'/f'{label}.ena.gff3.gz').exists())
        fields = dict(l.split('\t',1) for l in (self.root/'out/manifest.txt').read_text().splitlines())
        self.assertEqual(fields['NAME'], label+' decoupled annotation')
        self.assertIn('Model: Insecta',fields['DESCRIPTION'])
        report=json.loads((self.root/'out/preparation.json').read_text())
        self.assertEqual(report['submission_metadata']['predictor_model'], 'Insecta')
        self.assertEqual(report['workflow']['source_detection']['model'], 'Insecta')

    def test_omitted_model_has_no_empty_component(self):
        s = g.resolve_settings(self.settings)
        self.assertEqual(s['label'], 'GCA_123456789.1.ANNEVO.2.2.3')
        self.assertNotIn('Model:',s['description'])

    def test_custom_description_is_preserved(self):
        s=g.resolve_settings({**self.settings,'predictor_model':'Insecta','description':'Custom annotation.'})
        self.assertEqual(s['description'], 'Custom annotation.')

    def test_explicit_generic_override(self):
        s=g.resolve_settings({**self.settings,'source_format':'generic'})
        self.assertEqual(s['source_format'], 'generic')

    def test_length_comment_whitespace_variant(self):
        sizes=g.annevo_lengths(['# --- Prediction on sequence number 1 ( length=100 , name=seq1 ) ---'])
        self.assertEqual(sizes, {'seq1':100})

    def test_version_header_variants(self):
        self.source.write_text(self.source.read_text().replace('ANNEVO v2.2.3', 'ANNEVO version 2.2.3'))
        self.assertEqual(g.detect_source(self.source), ('ANNEVO','2.2.3'))

    def test_batch_model_field(self):
        tsv=self.root/'batch.tsv'
        tsv.write_text('input\tassembly_accession\tlocus_tag_prefix\tpredictor_model\ninput.gff\tGCA_123456789.1\tTEST\tInsecta\n')
        import contextlib,io
        with contextlib.redirect_stdout(io.StringIO()):
            code=g.batch_prepare(tsv, {'out_root':str(self.root/'batch')}, self.root/'summary.json')
        self.assertEqual(code,0)
        self.assertTrue((self.root/'batch/GCA_123456789.1.ANNEVO.2.2.3.Insecta/manifest.txt').exists())

    def test_model_requires_version_even_with_no_inference(self):
        self.source.write_text(self.source.read_text().replace('ANNEVO v2.2.3.', 'ANNEVO.'))
        with self.assertRaisesRegex(ValueError,'predictor-model requires'):
            g.resolve_settings({**self.settings,'no_inference':True,'predictor_model':'Insecta'})

    def test_annotator_version_model_description(self):
        settings = {**self.settings, "predictor_model": "Insecta"}
        result = g.resolve_settings(settings)
        self.assertEqual(result["description"], "Ab initio genome annotation prepared for ENA decoupled submission. Annotator: ANNEVO Version: 2.2.3 Model: Insecta")
        self.assertEqual(result["inference"], "ab initio prediction:ANNEVO:2.2.3")
        without_model = g.resolve_settings(self.settings)
        self.assertEqual(without_model["description"], "Ab initio genome annotation prepared for ENA decoupled submission. Annotator: ANNEVO Version: 2.2.3")
