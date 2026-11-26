import importlib.util
import tempfile
import unittest
from pathlib import Path
from argparse import Namespace

spec = importlib.util.spec_from_file_location('prep', Path(__file__).resolve().parents[1] / 'ena_gffprep.py')
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

BASE = '''##gff-version 3
##sequence-region seq1 1 1000
seq1\ttest\tgene\t1\t900\t.\t+\t.\tID=g1
seq1\ttest\tmRNA\t1\t900\t.\t+\t.\tID=t1;Parent=g1
seq1\ttest\tCDS\t10\t100\t.\t+\t0\tID=c1;Parent=t1;product=test%2C protein
seq1\ttest\tCDS\t201\t299\t.\t+\t2\tID=c1;Parent=t1;product=test%2C protein
'''

class Tests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.d = Path(self.tmp.name)
        self.n = 0
    def runprep(self, text=BASE, **kw):
        self.n += 1
        d = self.d / str(self.n)
        d.mkdir()
        (d/'in.gff').write_text(text)
        a = Namespace(input=str(d/'in.gff'), output=str(d/'out.gff.gz'),
            write_locus_tag_map=str(d/'tags.tsv'), report=str(d/'report.json'),
            locus_tag_prefix='TEST', width=6, start=1, step=1, taxon_id=9606,
            translation_table=1, locus_tag_map=None, assembly_report=None,
            accession_column='genbank', fai=None, strict_warnings=False,
            default_product=None, assembly_accession=None,bioproject='PRJEB123456', manifest=None,
            name='Test', description='Test')
        for k,v in kw.items():
            setattr(a,k,v)
        report = m.prepare(a)
        with m.opened(a.output) as f:
            out = f.read()
        return a, report, out
    def test_basic_and_escaping(self):
        a,r,o = self.runprep()
        self.assertEqual(o.count('locus_tag=TEST_000001'),4)
        self.assertIn('product=test%2C protein',o)
        self.assertEqual(r['genes'],1)
        self.assertEqual(o.count('transl_table=1'),2)
    def test_mapping_reserves_retired(self):
        a,_,_ = self.runprep()
        b,r,o = self.runprep(BASE.replace('g1','g2').replace('t1','t2').replace('c1','c2'), locus_tag_map=a.write_locus_tag_map)
        self.assertIn('TEST_000002',o)
        self.assertEqual(r['reserved_retired_gene_ids'],['g1'])
        _,_,o = self.runprep(locus_tag_map=b.write_locus_tag_map)
        self.assertIn('TEST_000001',o)
    def test_missing_parent(self):
        with self.assertRaisesRegex(m.Invalid,'Missing Parent'):
            self.runprep(BASE.replace('Parent=g1','Parent=absent'))
    def test_cycle(self):
        with self.assertRaises(m.Invalid):
            self.runprep(BASE.replace('ID=g1','ID=g1;Parent=t1'))
    def test_bounds(self):
        with self.assertRaisesRegex(m.Invalid,'beyond sequence'):
            self.runprep(BASE.replace('1 1000','1 100'))
    def test_partial(self):
        with self.assertRaisesRegex(m.Invalid,'partial'):
            self.runprep(BASE.replace('ID=g1','ID=g1;partial=true'))
    def test_shared_leaf(self):
        text=BASE + 'seq1\ttest\tmRNA\t1\t900\t.\t+\t.\tID=t2;Parent=g1\nseq1\ttest\texon\t1\t100\t.\t+\t.\tID=e1;Parent=t1,t2\n'
        _,r,o=self.runprep(text)
        self.assertEqual(r['shared_features_split'],1)
        self.assertIn('ID=e1.ena_parent1;Parent=t1',o)
        self.assertIn('ID=e1.ena_parent2;Parent=t2',o)
    def test_conflicting_tags(self):
        with self.assertRaisesRegex(m.Invalid,'Conflicting locus'):
            self.runprep(BASE.replace('ID=g1','ID=g1;locus_tag=TEST_000001').replace('ID=t1','ID=t1;locus_tag=TEST_000002'))
    def test_pseudo(self):
        _,_,o=self.runprep(BASE.replace('ID=g1','ID=g1;pseudo=true'))
        self.assertNotIn('transl_table=',o)
    def test_fasta_rejected(self):
        with self.assertRaisesRegex(m.Invalid,'Embedded FASTA'):
            self.runprep(BASE+'##FASTA\n>seq1\nACGT\n')
    def test_report_remap(self):
        p=self.d/'GCA_123456789.1_sequence_report.txt'
        p.write_text('accession\tsequence-name\tsequence-length\tsequence-role\treplicon-name\treplicon-type\tassembly-unit\nAB123456.1\tseq1\t1000\tassembled-molecule\t1\tChromosome\tPrimary Assembly\n')
        _,_,o=self.runprep(assembly_report=str(p),assembly_accession="GCA_123456789.1")
        self.assertIn('##sequence-region AB123456.1 1 1000',o)
        self.assertNotIn('\nseq1\t',o)
    def test_idempotent(self):
        a,_,o=self.runprep()
        _,_,o2=self.runprep(o,locus_tag_map=a.write_locus_tag_map)
        self.assertEqual(o,o2)
    def test_orphan(self):
        with self.assertRaisesRegex(m.Invalid,'Orphan'):
            self.runprep(BASE.replace('ID=t1;Parent=g1','ID=t1'))
    def test_manifest_preflight(self):
        with self.assertRaisesRegex(m.Invalid,'assembly-accession'):
            self.runprep(manifest=str(self.d/'manifest.txt'))

    def test_annevo_grouping(self):
        text=BASE.replace('ID=c1;','ID=t1.CDS.1;',1).replace('ID=c1;','ID=t1.CDS.2;',1)
        # Identical descriptive attributes across segments are retained.
        _,report,out=self.runprep(text,source_format='annevo')
        self.assertEqual(report['adapter']['cds_rows_renamed'],2)
        self.assertEqual(out.count('ID=t1.CDS;Parent=t1'),2)
        self.assertIn('original_id=t1.CDS.1',out)
        self.assertIn('original_id=t1.CDS.2',out)
    def test_existing_logical_id_is_source_independent(self):
        _,_,out=self.runprep(source_format='annevo')
        self.assertEqual(out.count('ID=c1;Parent=t1'),2)
    def test_annevo_comment_lengths(self):
        text=BASE.replace('##sequence-region seq1 1 1000','# ----- prediction on sequence number 1 (length = 1000, name = seq1) -----')
        text=text.replace('ID=c1;','ID=t1.CDS.1;',1).replace('ID=c1;','ID=t1.CDS.2;',1)
        _,report,out=self.runprep(text,source_format='annevo',use_annevo_lengths=True)
        self.assertIn('##sequence-region seq1 1 1000',out)
        self.assertEqual(report['sequence_lengths_source'],'annevo_comments')
    def test_annevo_length_disagreement(self):
        text=BASE+'# ----- prediction on sequence number 1 (length = 999, name = seq1) -----\n'
        text=text.replace('ID=c1;','ID=t1.CDS.1;',1).replace('ID=c1;','ID=t1.CDS.2;',1)
        with self.assertRaisesRegex(m.Invalid,'length mismatch'):
            self.runprep(text,source_format='annevo',use_annevo_lengths=True)
    def test_annevo_changed_product_rejected(self):
        text=BASE.replace('ID=c1;','ID=t1.CDS.1;',1).replace('ID=c1;','ID=t1.CDS.2;',1)
        text=text.replace('product=test%2C protein','product=other',1)
        with self.assertRaisesRegex(m.Invalid,'Conflicting product.*CDS attributes'):
            self.runprep(text,source_format='annevo')
    def test_audit_phase_checks(self):
        p=self.d/'audit.gff'
        # First segment is 91 bases: phase of following segment must be 2.
        p.write_text(BASE)
        report=m.audit(p)
        self.assertEqual(report['phase_continuity_errors'],[])
        p.write_text(BASE.replace('\t+\t2\t','\t+\t1\t'))
        report=m.audit(p)
        self.assertEqual(len(report['phase_continuity_errors']),1)
    def test_annevo_output_idempotent(self):
        text=BASE.replace('ID=c1;','ID=t1.CDS.1;',1).replace('ID=c1;','ID=t1.CDS.2;',1)
        a,_,out=self.runprep(text,source_format='annevo')
        _,_,out2=self.runprep(out,source_format='annevo',locus_tag_map=a.write_locus_tag_map)
        self.assertEqual(out,out2)

    def test_species_conflict(self):
        with self.assertRaisesRegex(m.Invalid,'species'):
            self.runprep('##species https://www.ncbi.nlm.nih.gov/datasets/taxonomy/706635/\n'+BASE)
    def test_species_match(self):
        _,_,out=self.runprep('##species https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/\n'+BASE)
        self.assertEqual(out.count('##species '),1)
    def test_tab_separated_directive(self):
        _,_,out=self.runprep(BASE.replace('##sequence-region seq1 1 1000','##sequence-region\tseq1\t1\t1000'))
        self.assertIn('##sequence-region seq1 1 1000',out)
    def test_reject_gff2(self):
        with self.assertRaisesRegex(m.Invalid,'only GFF3'):
            self.runprep(BASE.replace('##gff-version 3','##gff-version 2'))
    def test_empty_parent_list_element(self):
        with self.assertRaisesRegex(m.Invalid,'Empty list member'):
            self.runprep(BASE.replace('Parent=t1','Parent=t1,'))
    def test_bad_escape(self):
        with self.assertRaisesRegex(m.Invalid,'percent'):
            self.runprep(BASE.replace('test%2C protein','test%protein'))
    def test_reserved_escaping_roundtrip(self):
        text=BASE.replace('test%2C protein','x%25%3B%3D%26%2Cy')
        a,_,out=self.runprep(text)
        _,_,out2=self.runprep(out,locus_tag_map=a.write_locus_tag_map)
        self.assertEqual(out,out2)
    def test_whitespace_ID(self):
        with self.assertRaisesRegex(m.Invalid,'whitespace'):
            self.runprep(BASE.replace('ID=g1','ID=g 1'))
    def test_fai_missing_id_not_hidden_by_region(self):
        p=self.d/'wrong.fai';p.write_text('other\t1000\t0\t80\t81\n')
        with self.assertRaisesRegex(m.Invalid,'absent from external'):
            self.runprep(fai=str(p))
    def test_fai_duplicate_id(self):
        p=self.d/'duplicate.fai';p.write_text('seq1\t1000\nseq1\t1000\n')
        with self.assertRaisesRegex(m.Invalid,'Duplicate FAI'):
            self.runprep(fai=str(p))
    def test_fai_length_mismatch(self):
        p=self.d/'short.fai';p.write_text('seq1\t999\n')
        with self.assertRaisesRegex(m.Invalid,'length disagrees'):
            self.runprep(fai=str(p))
    def test_invalid_genetic_code(self):
        with self.assertRaisesRegex(m.Invalid,'Translation table'):
            self.runprep(translation_table=8)
    def test_invalid_existing_genetic_code(self):
        with self.assertRaisesRegex(m.Invalid,'Invalid transl_table'):
            self.runprep(BASE.replace('ID=c1;','transl_table=999;ID=c1;'))
    def test_conflicting_discontinuous_product(self):
        with self.assertRaisesRegex(m.Invalid,'Conflicting product'):
            self.runprep(BASE.replace('product=test%2C protein','product=other',1))
    def test_preserve_discontinuous_translation_code(self):
        text=BASE.replace('ID=c1;','transl_table=5;ID=c1;',1)
        _,_,out=self.runprep(text)
        self.assertEqual(out.count('transl_table=5'),2)
        self.assertNotIn('transl_table=1',out)
    def test_phase_review_in_preparation(self):
        text=BASE.replace('\t+\t2\t','\t+\t1\t')
        _,report,_=self.runprep(text)
        self.assertEqual(len(report['biological_review']['phase_continuity_errors']),1)
    def test_minus_phase_review(self):
        # On reverse strand the 201..299 interval is first; 99 nt -> next phase 0.
        text=BASE.replace('\t+\t','\t-\t').replace('\t-\t2\t','\t-\t0\t')
        _,report,_=self.runprep(text)
        self.assertEqual(report['biological_review']['phase_continuity_errors'],[])
    def test_missing_product_is_note_not_strict_error(self):
        text=BASE.replace('\t10\t100\t','\t10\t99\t').replace('\t+\t2\t','\t+\t0\t').replace(';product=test%2C protein','')
        _,report,out=self.runprep(text,strict_warnings=True)
        self.assertEqual(report['warnings'],[])
        self.assertTrue(any('lack product' in note for note in report['notes']))
        self.assertNotIn('product=',out)
    def test_output_subtree_contiguous(self):
        text=BASE+'seq1\ttest\tgene\t2\t950\t.\t+\t.\tID=g2\n'
        _,_,out=self.runprep(text)
        self.assertLess(out.rfind('ID=c1'),out.index('ID=g2'))
    def test_pseudogene_flags_propagated(self):
        _,_,out=self.runprep(BASE.replace('ID=g1','ID=g1;pseudo=true'))
        self.assertEqual(out.count('pseudo=true'),4)
    def test_transcript_id_reused(self):
        text=BASE.replace('ID=t1;','transcript_id=shared;ID=t1;')+'seq1\ttest\tmRNA\t1\t900\t.\t+\t.\tID=t2;Parent=g1;transcript_id=shared\n'
        with self.assertRaisesRegex(m.Invalid,'multiple logical features'):
            self.runprep(text)
    def test_special_interval_order_refused(self):
        with self.assertRaisesRegex(m.Invalid,'specialist handling'):
            self.runprep(BASE.replace('ID=c1;','part=1;ID=c1;'))
    def test_gap_required_length(self):
        text=BASE+'seq1\ttest\tgap\t950\t959\t.\t.\t.\tID=gap1;gap_type=within scaffold\n'
        with self.assertRaisesRegex(m.Invalid,'estimated_length'):
            self.runprep(text)
    def test_recoded_codon_phase_and_bounds(self):
        text=BASE+'seq1\ttest\trecoded_codon\t30\t32\t.\t+\t0\tID=codon;Parent=c1;codon_redefined=selenocysteine\n'
        _,_,out=self.runprep(text)
        self.assertIn('codon_redefined=selenocysteine',out)
        with self.assertRaisesRegex(m.Invalid,'intron'):
            self.runprep(text.replace('\t30\t32\t','\t150\t152\t'))
    def test_map_padding_collision(self):
        p=self.d/'map.tsv';p.write_text('gene_id\tlocus_tag\ng1\tTEST_1\ng2\tTEST_000001\n')
        with self.assertRaisesRegex(m.Invalid,'Ambiguous numeric'):
            self.runprep(locus_tag_map=str(p))
    def test_report_metadata_and_hashes(self):
        _,report,_=self.runprep()
        self.assertEqual(report['submission_metadata']['taxon_id'],9606)
        self.assertEqual(len(report['input_files_sha256']),1)
    def test_manifest_relative_and_tabs(self):
        dest=self.d/'manifest.txt'
        a,_,_=self.runprep(manifest=str(dest),assembly_accession='GCA_948150665.1')
        fields=dict(line.split('\t',1) for line in dest.read_text().splitlines())
        self.assertEqual(fields['ANALYSIS_TYPE'],'DECOUPLED_ANNOTATION')
        self.assertEqual(fields['PRIMARY_ID'],'GCA_948150665.1')
        self.assertEqual((dest.parent/fields['GFF3']).resolve(),Path(a.output).resolve())
    def test_manifest_bioproject_study(self):
        dest=self.d/'study_manifest.txt'
        a,_,_=self.runprep(manifest=str(dest),assembly_accession='GCA_948150665.1',bioproject='PRJEB123456')
        fields=dict(line.split('\t',1) for line in dest.read_text().splitlines())
        self.assertEqual(fields['STUDY'],'PRJEB123456')
        self.assertEqual(fields['PRIMARY_ID'],'GCA_948150665.1')
        self.assertEqual((dest.parent/fields['GFF3']).resolve(),Path(a.output).resolve())

    def test_commit_failure_rolls_back(self):
        from unittest.mock import patch
        real_link=m.os.link
        calls=[]
        def failing_link(source,dest):
            calls.append(dest)
            if len(calls)==2:
                raise OSError('Simulated output publication failure')
            return real_link(source,dest)
        with patch.object(m.os,'link',side_effect=failing_link):
            with self.assertRaises(OSError):
                self.runprep()
        self.assertFalse((self.d/'1'/'out.gff.gz').exists())
        self.assertFalse((self.d/'1'/'tags.tsv').exists())
        self.assertEqual(list((self.d/'1').glob('.ena_gffprep-*')),[])

if __name__ == '__main__':
    unittest.main()
