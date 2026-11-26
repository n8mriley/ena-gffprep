import unittest
import ena_ontology as so
import test_gffprep as base


class OntologyTests(unittest.TestCase):
    setUp = base.Tests.setUp
    runprep = base.Tests.runprep

    def test_pinned_graph(self):
        self.assertEqual(so.DATA['version'], '2026-08-07')
        self.assertEqual(len(so.TERMS), 2620)
        self.assertTrue(so.is_a('mt_gene', 'gene'))
        self.assertTrue(so.is_a('interior_coding_exon', 'exon'))
        self.assertEqual(so.resolve('SO:0000358'), 'SO:0000104')

    def test_nonfeature(self):
        with self.assertRaisesRegex(base.m.Invalid, 'not an SO sequence feature'):
            self.runprep(base.BASE.replace('\tgene\t','\tprotein_coding\t'))

    def test_obsolete(self):
        with self.assertRaisesRegex(base.m.Invalid, 'obsolete SO term.*SO:0001556'):
            self.runprep(base.BASE.replace('\tgene\t','\tSO:0000053\t'))

    def test_valid_but_unmapped(self):
        with self.assertRaisesRegex(base.m.Invalid, 'local INSDC vocabulary.*sequence_feature'):
            self.runprep(base.BASE.replace('\tgene\t','\tmt_gene\t'))

    def test_relationship_evidence(self):
        self.assertEqual(so.relationship('riboswitch','mRNA'), 'ontology_supported')
        self.assertEqual(so.relationship('gene','exon'), 'not_established')
        self.assertEqual(so.relationship('unknown','gene'), 'unresolved_term')

    def test_report(self):
        _, report, _ = self.runprep(base.BASE)
        self.assertEqual(report['ontology']['version'], so.DATA['version'])
        self.assertEqual(report['ontology']['gene_subtypes'], ['gene'])
        self.assertGreater(sum(report['ontology']['parent_relationships_by_logical_feature'].values()), 0)
