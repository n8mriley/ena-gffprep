import unittest
from ena_gffprep import parser, resolve_settings

class BioProjectTests(unittest.TestCase):
    def test_cli_flag(self):
        self.assertEqual(parser().parse_args(['--bioproject','PRJEB123456']).bioproject, 'PRJEB123456')

    def test_setting_recognized_in_config_and_batch(self):
        from ena_gffprep import SETTING_KEYS
        self.assertIn('bioproject', SETTING_KEYS)


class MandatoryStudyTests(unittest.TestCase):
    def setUp(self):
        import tempfile
        from pathlib import Path
        from test_gffprep import BASE
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        source = self.root / 'input.gff'
        source.write_text(BASE)
        self.settings = dict(input=str(source), assembly_accession='GCA_123456789.1',
                             locus_tag_prefix='TEST', no_inference=True,
                             output_dir=str(self.root / 'out'))

    def test_missing_or_invalid_study(self):
        for value in (None, '', 'PRJEB', 'ERP123', 'PRJEB12\n', 123, 'PRJEB１２'):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'bioproject'):
                resolve_settings(dict(self.settings, bioproject=value))
        self.assertFalse((self.root / 'out').exists())

    def test_valid_prefixes(self):
        for value in ('PRJEB123', 'PRJNA123', 'PRJDB123'):
            self.assertEqual(resolve_settings(dict(self.settings, bioproject=value))['bioproject'], value)

    def test_config_study_and_manifest(self):
        import json
        from ena_gffprep import configuration, prepare_package
        config = self.root / 'config.json'
        config.write_text(json.dumps({'bioproject': 'PRJEB123'}))
        prepare_package(dict(self.settings, **configuration(config)))
        self.assertIn('STUDY\tPRJEB123\n', (self.root / 'out/manifest.txt').read_text())
        report = json.loads((self.root / 'out/preparation.json').read_text())
        self.assertEqual(report['submission_metadata']['bioproject'], 'PRJEB123')

    def test_batch_missing_and_row_study(self):
        import contextlib, io, json
        from ena_gffprep import batch_prepare
        batch = self.root / 'jobs.tsv'
        batch.write_text('input\tbioproject\tlocus_tag_prefix\ninput.gff\t\tTESA\ninput.gff\tPRJEB123\tTESB\n')
        defaults = {k:v for k,v in self.settings.items() if k not in {'input','output_dir','locus_tag_prefix'}}
        defaults['out_root'] = str(self.root / 'batch')
        summary = self.root / 'summary.json'
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(batch_prepare(batch, defaults, summary), 1)
        rows = json.loads(summary.read_text())['rows']
        self.assertEqual([r['status'] for r in rows], ['failed','prepared'])
        self.assertIn('bioproject', rows[0]['error'])

    def test_audit_exempt(self):
        from ena_gffprep import main
        report = self.root / 'audit.json'
        self.assertEqual(main(['--input', self.settings['input'], '--audit-only', '--report', str(report)]), 0)
        self.assertTrue(report.exists())

    def test_direct_api_requires_study(self):
        from test_gffprep import Tests, m
        fixture = Tests()
        fixture.setUp()
        try:
            with self.assertRaisesRegex(m.Invalid, 'bioproject'):
                fixture.runprep(bioproject=None)
        finally:
            fixture.doCleanups()

if __name__ == '__main__':
    unittest.main()
