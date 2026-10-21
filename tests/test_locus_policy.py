"""Regression checks for DDBJ locus-tag guidance and stable ENA allocation."""
import unittest
import test_gffprep as fixtures
BASE, m = fixtures.BASE, fixtures.m

class LocusPolicy(unittest.TestCase):
    setUp = fixtures.Tests.setUp
    runprep = fixtures.Tests.runprep

    def mapfile(self, rows):
        path = self.d / 'mapping.tsv'
        path.write_text('gene_id\tlocus_tag\n' + rows)
        return str(path)

    def test_alphanumeric_tags_preserved_and_roundtrip(self):
        a, _, out = self.runprep(locus_tag_map=self.mapfile('g1\tTEST_A0001\nold\tTEST_r1112\n'))
        self.assertIn('locus_tag=TEST_A0001', out)
        self.assertIn('old\tTEST_r1112\treserved', m.Path(a.write_locus_tag_map).read_text())
        _, _, again = self.runprep(out, locus_tag_map=a.write_locus_tag_map)
        self.assertEqual(out, again)

    def test_mixed_reserved_suffixes_high_water(self):
        _, _, out = self.runprep(locus_tag_map=self.mapfile('retired\tTEST_004321\nother\tTEST_Z9999\n'))
        self.assertIn('locus_tag=TEST_004322', out)

    def test_decimal_version_rejected(self):
        with self.assertRaisesRegex(m.Invalid, 'no decimal versions'):
            self.runprep(locus_tag_map=self.mapfile('g1\tTEST_0020.1\n'))

    def test_wrong_prefix_rejected(self):
        with self.assertRaisesRegex(m.Invalid, 'registered prefix'):
            self.runprep(locus_tag_map=self.mapfile('g1\tOTHER_A1\n'))

    def test_alphanumeric_collision_rejected(self):
        with self.assertRaisesRegex(m.Invalid, 'reused by different genes'):
            self.runprep(locus_tag_map=self.mapfile('g1\tTEST_A1\ng2\tTEST_A1\n'))

    def test_repeat_region_not_tagged(self):
        text = BASE + 'seq1\ttest\trepeat_region\t20\t40\t.\t+\t.\tID=r1;Parent=g1\n'
        _, _, out = self.runprep(text)
        repeat = next(line for line in out.splitlines() if '\trepeat_region\t' in line)
        self.assertNotIn('locus_tag=', repeat)

    def test_existing_repeat_tag_rejected(self):
        text = BASE + 'seq1\ttest\trepeat_region\t20\t40\t.\t+\t.\tID=r1;Parent=g1;locus_tag=TEST_000001\n'
        with self.assertRaisesRegex(m.Invalid, 'must not carry locus_tag'):
            self.runprep(text)

    def test_numbering_across_sequences_and_gene_types(self):
        text = BASE.replace('##sequence-region seq1 1 1000', '##sequence-region seq1 1 1000\n##sequence-region seq2 1 1000')
        text += 'seq2\ttest\tgene\t1\t80\t.\t+\t.\tID=g2\nseq2\ttest\ttRNA\t1\t80\t.\t+\t.\tID=r2;Parent=g2\n'
        _, _, out = self.runprep(text, width=5)
        self.assertIn('ID=g1;locus_tag=TEST_00001', out)
        self.assertIn('ID=g2;locus_tag=TEST_00002', out)
        self.assertIn('ID=r2;Parent=g2;locus_tag=TEST_00002', out)

    def test_digit_in_prefix_supported(self):
        _, _, out = self.runprep(locus_tag_prefix='A1C', width=5)
        self.assertIn('locus_tag=A1C_00001', out)

    def test_rna_components_share_gene_tag(self):
        for kind in ('tRNA', 'rRNA', 'ncRNA', 'tmRNA', 'precursor_RNA', 'prim_transcript', 'misc_RNA'):
            with self.subTest(kind=kind):
                text = '##gff-version 3\n##sequence-region seq1 1 1000\n'
                text += 'seq1\ttest\tgene\t1\t80\t.\t+\t.\tID=g1\n'
                text += f'seq1\ttest\t{kind}\t1\t80\t.\t+\t.\tID=r1;Parent=g1\n'
                _, _, out = self.runprep(text)
                self.assertIn('ID=r1;Parent=g1;locus_tag=TEST_000001', out)
