"""Submission-profile regressions for the initial conformance update."""
import unittest
import test_gffprep as base


class Conformance(unittest.TestCase):
    setUp = base.Tests.setUp
    runprep = base.Tests.runprep

    def test_species_host(self):
        with self.assertRaisesRegex(base.m.Invalid, 'species URI'):
            self.runprep('##species https://example.org/taxonomy/9606/\n'+base.BASE)

    def test_primary_pragma(self):
        _, report, output = self.runprep('#!transl_table primary:6\n'+base.BASE)
        self.assertEqual(output.count('transl_table=6'), 2)
        self.assertEqual(report['translation_table_sources_by_CDS_row'], {'primary_pragma': 2})

    def test_cds_override(self):
        _, _, output = self.runprep('#!transl_table primary:6\n'+base.BASE.replace('ID=c1;', 'ID=c1;transl_table=5;'))
        self.assertEqual(output.count('transl_table=5'), 2)

    def test_compartment_pragma_refused(self):
        with self.assertRaisesRegex(base.m.Invalid, 'compartment'):
            self.runprep('#!transl_table primary:1,mitochondrion:2\n'+base.BASE)

    def test_internal_partial_both_strands(self):
        for strand, flag in [('+','end'), ('-','start')]:
            text = base.BASE.replace('\t+\t', '\t'+strand+'\t').replace('ID=c1;', 'ID=c1;partial='+flag+';', 1)
            with self.assertRaisesRegex(base.m.Invalid, 'Internal partial'):
                self.runprep(text)

    def test_terminal_partial(self):
        self.runprep(base.BASE.replace('ID=c1;', 'ID=c1;partial=start;', 1))

    def test_amino_acid(self):
        with self.assertRaisesRegex(base.m.Invalid, 'amino acid'):
            self.runprep(base.BASE+'seq1\ttest\trecoded_codon\t30\t32\t.\t+\t.\tID=codon;Parent=c1;codon_redefined=banana\n')

    def test_dbxref_syntax(self):
        with self.assertRaisesRegex(base.m.Invalid, 'Dbxref'):
            self.runprep(base.BASE.replace('ID=g1','ID=g1;Dbxref=missingcolon'))

    def test_unknown_directive(self):
        _, report, output = self.runprep('##unknown-directive x\n'+base.BASE)
        self.assertNotIn('##unknown-directive', output)
        self.assertIn('##unknown-directive x', report['removed_comments'])

    def test_exceptional_stop(self):
        row = 'seq1\ttest\tstop_codon\t30\t32\t.\t+\t.\tID=stop;Parent=c1\n'
        with self.assertRaisesRegex(base.m.Invalid, 'biological CDS end'):
            self.runprep(base.BASE+row)
        self.runprep(base.BASE+row.replace('30\t32','298\t299'))

    def test_unknown_feature(self):
        with self.assertRaisesRegex(base.m.Invalid, 'local INSDC vocabulary'):
            self.runprep(base.BASE.replace('\tgene\t', '\tfictional_feature\t'))

    def test_so_accession(self):
        _, _, output = self.runprep(base.BASE.replace('\tCDS\t', '\tSO:0000316\t'))
        self.assertEqual(output.count('\tCDS\t'), 2)

    def test_so_accession_phase(self):
        with self.assertRaisesRegex(base.m.Invalid, 'invalid phase'):
            self.runprep(base.BASE.replace('\tCDS\t', '\tSO:0000316\t').replace('\t0\tID=c1', '\t.\tID=c1'))

    def test_dbxref_registry(self):
        self.runprep(base.BASE.replace('ID=g1', 'ID=g1;Dbxref=GeneID:123,GO:GO:123,ATCC(in host):123'))
        for tag in ('FakeDB', 'geneid'):
            with self.assertRaisesRegex(base.m.Invalid, 'unrecognized Dbxref'):
                self.runprep(base.BASE.replace('ID=g1', 'ID=g1;Dbxref='+tag+':123'))

    def test_pseudogene_subtype(self):
        self.runprep(base.BASE.replace('ID=g1', 'ID=g1;pseudogene=processed'))
        with self.assertRaisesRegex(base.m.Invalid, 'pseudogene subtype'):
            self.runprep(base.BASE.replace('ID=g1', 'ID=g1;pseudogene=fictional'))

    def test_specific_pseudogene_feature(self):
        _, report, output = self.runprep(base.BASE.replace('\tgene\t', '\tprocessed_pseudogene\t'))
        self.assertEqual(report['genes'], 1)
        self.assertIn('pseudogene=processed', output)
        self.assertNotIn('transl_table=', output)

    def test_invalid_cds_parent(self):
        text = base.BASE.replace('ID=c1;Parent=t1', 'ID=c1;Parent=e1')
        text += 'seq1\ttest\texon\t1\t900\t.\t+\t.\tID=e1;Parent=t1\n'
        with self.assertRaisesRegex(base.m.Invalid, 'Invalid parent type'):
            self.runprep(text)

    def test_compartment_codes(self):
        acc = 'GCA_123456789.1'
        path = self.d / (acc+'_sequence_report.txt')
        path.write_text('accession\tsequence-name\tsequence-length\tsequence-role\treplicon-name\treplicon-type\tassembly-unit\n'
                        'OX123456.1\tseq1\t1000\tassembled-molecule\t1\tChromosome\tPrimary Assembly\n'
                        'OX123457.1\tseq2\t1000\tassembled-molecule\tMT\tMitochondrion\tnon-nuclear\n')
        second = '\n'.join(line.replace('seq1','seq2').replace('g1','g2').replace('t1','t2').replace('c1','c2') for line in base.BASE.splitlines()[1:])+'\n'
        text = '#!transl_table primary:1,mitochondrion:2\n'+base.BASE+second
        _, report, output = self.runprep(text, assembly_accession=acc, sequence_metadata=str(path))
        self.assertEqual(report['translation_table_sources_by_CDS_row'], {'primary_pragma':2,'mitochondrion_pragma':2})
        self.assertEqual(output.count('transl_table=2'),2)
        with self.assertRaisesRegex(base.m.Invalid, 'Missing or uncertain'):
            self.runprep(text.replace(',mitochondrion:2',''), assembly_accession=acc, sequence_metadata=str(path))

    def test_duplicate_compartment(self):
        with self.assertRaisesRegex(base.m.Invalid, 'Duplicate transl_table compartment'):
            self.runprep('#!transl_table primary:1,primary:2\n'+base.BASE)

    def test_compartment_explicit_override(self):
        _, report, output = self.runprep('#!transl_table primary:1,mitochondrion:2\n'+base.BASE.replace('ID=c1;', 'ID=c1;transl_table=5;'))
        self.assertEqual(report['translation_table_sources_by_CDS_row'], {'CDS_attribute':2})

    def test_slippage_overlap_preserved(self):
        text = base.BASE.replace('201\t299', '100\t299').replace('ID=c1;', 'ID=c1;exception=ribosomal_slippage;')
        _, report, output = self.runprep(text)
        self.assertIn('\t100\t299\t.\t+\t2\t', output)
        self.assertIn('c1', report['biological_review']['phase_check_skipped_exception_or_pseudogene'])
        _, _, minus = self.runprep(text.replace('\t+\t', '\t-\t'))
        self.assertIn('\t100\t299\t.\t-\t2\t', minus)
        with self.assertRaisesRegex(base.m.Invalid, 'Overlapping'):
            self.runprep(text.replace('exception=ribosomal_slippage;', ''))
        with self.assertRaisesRegex(base.m.Invalid, 'Overlapping'):
            self.runprep(text.replace('100\t299', '99\t299'))

    def test_split_recoded_codon(self):
        text = base.BASE + 'seq1\ttest\trecoded_codon\t99\t100\t.\t+\t.\tID=rc;Parent=c1;codon_redefined=selenocysteine\nseq1\ttest\trecoded_codon\t201\t201\t.\t+\t.\tID=rc;Parent=c1;codon_redefined=selenocysteine\n'
        _, _, output = self.runprep(text)
        self.assertEqual(output.count('codon_redefined=selenocysteine'), 2)
        with self.assertRaisesRegex(base.m.Invalid, 'Conflicting codon_redefined'):
            self.runprep(text.rsplit('selenocysteine',1)[0]+'methionine\n')

    def test_additional_rna_ownership(self):
        row = 'seq1\ttest\tSRP_RNA\t400\t500\t.\t+\t.\tID=rna;Parent=g1\n'
        _, _, output = self.runprep(base.BASE+row)
        self.assertIn('ID=rna;Parent=g1;locus_tag=', output)
        with self.assertRaisesRegex(base.m.Invalid, 'Orphan gene component'):
            self.runprep(base.BASE+row.replace(';Parent=g1',''))

    def test_pseudogenic_cds(self):
        text = base.BASE.replace('\tCDS\t', '\tpseudogenic_CDS\t')
        _, _, output = self.runprep(text)
        self.assertEqual(output.count('\tCDS\t'), 2)
        self.assertEqual(output.count('pseudo=true'), 2)
        self.assertNotIn('transl_table=', output)
        with self.assertRaisesRegex(base.m.Invalid, 'invalid phase'):
            self.runprep(text.replace('\t0\tID=c1', '\t.\tID=c1'))
