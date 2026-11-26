"""Exact user-supplied organelle exclusions, without remote classification."""
import unittest
from unittest.mock import patch
import ena_organelles as o
import ena_gffprep as cli
import test_gffprep as base
import test_organelles as fixtures

class ManualOrganelleTests(unittest.TestCase):
    setUp = base.Tests.setUp
    runprep = base.Tests.runprep

    def test_cli_modes(self):
        p = cli.parser()
        self.assertIsNone(p.parse_args([]).remove_organelles)
        self.assertIs(p.parse_args(['--remove-organelles']).remove_organelles, True)
        self.assertEqual(p.parse_args(['--remove-organelles','mt,pt']).remove_organelles, 'mt,pt')
        self.assertIs(p.parse_args(['--remove-organelles','--no-inference']).remove_organelles, True)

    def test_invalid_lists(self):
        for value in ['', 'mt,', ',mt', 'mt,,pt', 'mt,mt', 'mt name', 1, None]:
            with self.subTest(value=value), self.assertRaises(ValueError):
                o.manual_sequence_ids(value)
        self.assertEqual(o.manual_sequence_ids(' mt , pt '), {'mt','pt'})

    def test_manual_filter_without_network(self):
        with patch.object(o, 'load_metadata', side_effect=AssertionError('No lookup')):
            a, report, out = self.runprep(fixtures.TEXT, remove_organelles='mt', no_ena_lengths=True)
        self.assertNotIn('\nmt\t', out)
        self.assertNotIn('##sequence-region mt ', out)
        info = report['organelle_filter']
        self.assertEqual(info['source'], 'manual')
        self.assertEqual(info['removed_feature_rows'], 3)
        self.assertEqual(info['removed_genes'], 1)
        self.assertEqual(info['excluded_sequences'][0]['ena_sequence'], '')
        self.assertEqual(info['requested_sequence_ids'], ['mt'])

    def test_header_only_and_multiple(self):
        text = fixtures.TEXT + '##sequence-region pt 1 1000\n'
        _, report, out = self.runprep(text, remove_organelles='mt,pt', no_ena_lengths=True)
        self.assertNotIn('##sequence-region pt ', out)
        self.assertEqual(len(report['organelle_filter']['excluded_sequences']), 2)

    def test_absent_exact_version(self):
        with self.assertRaisesRegex(ValueError, 'absent from input'):
            self.runprep(fixtures.TEXT, remove_organelles='mt.1')

    def test_all_features_removed(self):
        with self.assertRaisesRegex(ValueError, 'remove all'):
            self.runprep(fixtures.TEXT, remove_organelles='mt,seq1')

    def test_cross_boundary_reference(self):
        text = fixtures.TEXT.replace('ID=t1;Parent=g1', 'ID=t1;Parent=gm')
        # Fixture IDs can vary; inject an explicit cross-boundary reference.
        text += 'seq1\ttest\texon\t1\t3\t.\t+\t.\tID=cross;Parent=gm\n'
        with self.assertRaisesRegex(ValueError, 'references an excluded'):
            self.runprep(text, remove_organelles='mt')

    def test_manual_rejects_organelle_lookup_options(self):
        with self.assertRaisesRegex(ValueError, 'cannot be combined'):
            self.runprep(fixtures.TEXT, remove_organelles='mt', organelle_metadata='unused.json')
