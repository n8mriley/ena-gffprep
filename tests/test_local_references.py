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
