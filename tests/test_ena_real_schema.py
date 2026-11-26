import unittest
from ena_organelles import parse_report, report_url

class ENARealSchema(unittest.TestCase):
    def test_observed_ena_rows(self):
        raw = ('accession\tsequence-name\tsequence-length\tsequence-role\treplicon-name\treplicon-type\tassembly-unit\n'
               'OX375756.1\t1\t70384342\tassembled-molecule\t1\tChromosome\tPrimary Assembly\n'
               'CANAHR010000001.1\tSUPER_X_unloc_1\t444687\tunlocalised-scaffold\tX\tChromosome\tPrimary Assembly\n'
               'OX375761.1\tMT\t18666\tassembled-molecule\tMT\tMitochondrion\tnon-nuclear\n').encode()
        rows = parse_report(raw,'GCA_947359355.1')
        self.assertEqual([r['length'] for r in rows],[70384342,444687,18666])
        self.assertEqual([r['organelle'] for r in rows],[False,False,True])
        self.assertEqual(report_url('GCA_947359355.1'),'https://ftp.ebi.ac.uk/pub/databases/ena/assembly/GCA_947/GCA_947359/GCA_947359355.1_sequence_report.txt')
    def test_unplaced_na(self):
        raw = ('accession\tsequence-name\tsequence-length\tsequence-role\treplicon-name\treplicon-type\tassembly-unit\n'
               'CBDEYJ010000002.1\tscaffold_20\t3620000\tunplaced-scaffold\tNA\tNA\tPrimary Assembly\n').encode()
        row = parse_report(raw,'GCA_965244925.1')[0]
        self.assertFalse(row['organelle']); self.assertTrue(row['unclassified'])
