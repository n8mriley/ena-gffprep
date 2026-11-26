import unittest
from unittest.mock import patch
from pathlib import Path
import ena_vocabularies as vocab
import test_gffprep as base

class LocalReferenceTests(unittest.TestCase):
    setUp = base.Tests.setUp
    runprep = base.Tests.runprep

    def test_missing_vocabulary(self):
        vocab.reference_data.cache_clear()
        try:
            with patch.object(Path, 'is_file', return_value=False):
                with self.assertRaisesRegex(ValueError, 'run tools/setup_references.py'):
                    vocab.reference_data()
        finally:
            vocab.reference_data.cache_clear()

    def test_reference_provenance(self):
        _, report, _ = self.runprep(base.BASE)
        self.assertEqual(report['local_reference_provenance']['jar_sha256'], vocab.reference_data()['jar_sha256'])
        self.assertEqual(report['local_reference_provenance']['dbxref_source_sha256'], vocab.reference_data()['dbxref_source_sha256'])

    def test_dbxref_compatibility_policy(self):
        raw = {'Bold', 'InterPro', 'Phytozome', 'Phytozme', 'TAIR'}
        self.assertEqual(vocab.effective_dbxref_tags(raw),
                         {'BOLD', 'Bold', 'InterPro', 'Interpro', 'Phytozome', 'TAIR'})
        self.assertIn('Phytozme', raw)
        self.assertEqual(vocab.effective_dbxref_tags({'Demo'}), {'Demo'})

    def test_dbxref_policy_in_preparation(self):
        for tag in ('BOLD', 'Bold', 'InterPro', 'Interpro', 'Phytozome', 'TAIR'):
            with self.subTest(tag=tag):
                self.runprep(base.BASE.replace('ID=g1', 'ID=g1;Dbxref='+tag+':123'))
        for tag in ('Phytozme', 'bold', 'interpro', 'NonexistentDatabase'):
            with self.subTest(tag=tag):
                with self.assertRaisesRegex(base.m.Invalid, 'unrecognized Dbxref'):
                    self.runprep(base.BASE.replace('ID=g1', 'ID=g1;Dbxref='+tag+':123'))
