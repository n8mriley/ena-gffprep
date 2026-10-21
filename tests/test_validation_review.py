import importlib.util
import tempfile
import unittest
from pathlib import Path
s=importlib.util.spec_from_file_location('review',Path(__file__).resolve().parents[1]/'ena_validation_review.py')
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
class ReviewTests(unittest.TestCase):
    def test_split_minus_cds_and_introns(self):
        with tempfile.TemporaryDirectory() as d:
            g=Path(d)/'a.gff';l=Path(d)/'report'
            g.write_text('##gff-version 3\ns\tx\tgene\t1\t100\t.\t-\t.\tID=g;locus_tag=TEST_1\ns\tx\tmRNA\t1\t100\t.\t-\t.\tID=t;Parent=g;locus_tag=TEST_1\ns\tx\tCDS\t1\t30\t.\t-\t0\tID=c;Parent=t;locus_tag=TEST_1\ns\tx\tCDS\t35\t64\t.\t-\t0\tID=c;Parent=t;locus_tag=TEST_1\n')
            l.write_text('ERROR: Violation of rule CDS_LENGTH on line 5: accession "s"\nWARNING: Violation of rule EXON_LENGTH: Violation of rule on line 4: test\n')
            summary,findings,short,introns=m.review(g,l)
            self.assertEqual(short[0]['coding_nt'],60);self.assertEqual(short[0]['segments'],2)
            self.assertEqual(findings[0]['gene_id'],'g');self.assertTrue(short[0]['reported_CDS_LENGTH'])
            self.assertEqual(introns[0]['intron_nt'],4)
            l.write_text('ERROR: Violation of rule CDS_LENGTH on line 99: test\n')
            with self.assertRaisesRegex(ValueError,'absent'):m.review(g,l)
            l.write_text('ERROR: Violation of rule CDS_LENGTH on line 5: accession "wrong"\n')
            with self.assertRaisesRegex(ValueError,'mismatch'):m.review(g,l)
