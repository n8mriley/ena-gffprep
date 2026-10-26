"""Source-independent rules, ambiguity guards and policy overrides."""
import gzip
import json
import tempfile
import unittest
from pathlib import Path
import ena_gffprep as g
import test_formats as format_fixtures
from test_organelles import ACC, record, MT

class StructuralTests(unittest.TestCase):
    fixture=format_fixtures.FormatTests.fixture
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);self.source=self.root/'input.gff';self.number=0
    def prep(self,text=None,**settings):
        self.number+=1;self.source.write_text('##gff-version 3\n##sequence-region s 1 1000\n'+(text if text is not None else self.fixture('UnknownPredictor')))
        config=dict(input=str(self.source),assembly_accession=ACC,locus_tag_prefix='TEST',no_inference=True,output_dir=str(self.root/str(self.number)),no_ncbi_lengths=True)
        config.update(settings);result=g.prepare_package(config);dest=Path(result['output_directory'])
        return gzip.open(next(dest.glob('*.gz')),'rt').read(),json.loads((dest/'preparation.json').read_text())
    def test_unknown_source_normalized_without_hint(self):
        out,r=self.prep()
        self.assertIn('\tmRNA\t',out);self.assertEqual(out.count('ID=t.CDS;'),2)
        self.assertEqual(r['source_format'],'generic');self.assertEqual(r['normalization']['engine'],'structural')
    def test_source_column_does_not_change_normalization(self):
        baseline=None;rules=None
        for source in ('UnknownPredictor','Helixer','Tiberius','Carbon-A','OrionGeno','AUGUSTUS','ANNEVO'):
            out,r=self.prep(self.fixture(source))
            neutral=out.replace('\t'+source+'\t','\tSOURCE\t')
            if baseline is None:baseline=neutral;rules=r['normalization']['rules']
            self.assertEqual(neutral,baseline);self.assertEqual(r['normalization']['rules'],rules)
    def test_named_compatibility_hints_do_not_select_engine(self):
        baseline=None
        for hint in ('auto','annevo','helixer','tiberius','braker','oriongeno','carbon-a','gff3'):
            out,r=self.prep(source_format=hint)
            if baseline is None:baseline=out
            self.assertEqual(out,baseline)
    def test_strict_mode_does_not_rewrite_types_or_cds_ids(self):
        out,r=self.prep(normalization='strict')
        self.assertIn('\ttranscript\t',out);self.assertIn('ID=a;Parent=t',out);self.assertNotIn('ID=t.CDS;',out)
    def test_explicit_generic_standard_override(self):
        out,r=self.prep(source_format='generic',normalization='standard')
        self.assertIn('\tmRNA\t',out);self.assertIn('ID=t.CDS;',out)
    def test_grouping_preserve_keeps_segment_ids(self):
        out,r=self.prep(cds_grouping='preserve')
        self.assertIn('\tmRNA\t',out);self.assertIn('ID=a;Parent=t',out);self.assertNotIn('ID=t.CDS;',out)
    def test_numbered_singleton_normalizes_for_every_source(self):
        text=self.fixture('Unknown').splitlines();text=text[:3];text[2]=text[2].replace('ID=a;','ID=t.CDS.1;')
        out,r=self.prep('\n'.join(text)+'\n')
        self.assertIn('ID=t.CDS;',out);self.assertIn('original_id=t.CDS.1',out)
    def test_existing_logical_group_retained(self):
        out,r=self.prep(self.fixture().replace('ID=a;','ID=logical;').replace('ID=b;','ID=logical;'))
        self.assertEqual(out.count('ID=logical;'),2)
    def test_distinct_proteins_preserved_as_products(self):
        text=self.fixture().replace('ID=a;Parent=t','ID=a;Parent=t;protein_id=P1').replace('ID=b;Parent=t','ID=b;Parent=t;protein_id=P2')
        out,r=self.prep(text)
        self.assertIn('ID=a;',out);self.assertIn('ID=b;',out)
        self.assertEqual(r['normalization']['rules']['explicit_product_groups_preserved'],2)
    def test_single_product_cannot_override_distinct_proteins(self):
        text=self.fixture().replace('ID=a;Parent=t','ID=a;Parent=t;protein_id=P1').replace('ID=b;Parent=t','ID=b;Parent=t;protein_id=P2')
        with self.assertRaisesRegex(ValueError,'Conflicting protein_id'):self.prep(text,cds_grouping='single-product')
    def test_multiple_multipart_groups_need_explicit_policy(self):
        text=self.fixture()+'s\tTiberius\tCDS\t40\t48\t.\t+\t0\tID=a;Parent=t\n'
        with self.assertRaisesRegex(ValueError,'Ambiguous existing multipart'):self.prep(text)
        out,r=self.prep(text,cds_grouping='preserve');self.assertIn('ID=a;',out)
        out,r=self.prep(text,cds_grouping='single-product');self.assertEqual(out.count('ID=t.CDS;'),3)
    def test_referenced_segments_are_not_renamed(self):
        with self.assertRaisesRegex(ValueError,'Referenced CDS segment'):
            self.prep(self.fixture()+'s\tTiberius\tstop_codon\t7\t9\t.\t+\t.\tID=stop;Parent=a\n')
    def test_unknown_metadata_conflict_not_silently_merged(self):
        text=self.fixture().replace('ID=a;Parent=t','ID=a;Parent=t;custom=A').replace('ID=b;Parent=t','ID=b;Parent=t;custom=B')
        with self.assertRaisesRegex(ValueError,'per-segment metadata'):self.prep(text)
    def test_segment_cds_type_preserved(self):
        text=self.fixture().replace('ID=a;Parent=t','ID=a;Parent=t;cds_type=initial').replace('ID=b;Parent=t','ID=b;Parent=t;cds_type=terminal')
        out,r=self.prep(text);self.assertIn('cds_type=initial',out);self.assertIn('cds_type=terminal',out)
    def test_existing_original_identifiers_retained(self):
        text=self.fixture().replace('ID=a;Parent=t','ID=a;Parent=t;original_id=older')
        out,r=self.prep(text);self.assertIn('original_id=older,a',out)
    def test_anonymous_cds_get_shared_logical_id(self):
        text=self.fixture().replace('ID=a;','').replace('ID=b;','')
        out,r=self.prep(text);self.assertEqual(out.count('ID=t.CDS'),2)
        self.assertEqual(r['normalization']['rules']['cds_rows_assigned_ID'],2)
    def test_mixed_anonymous_identified_segments_refused(self):
        with self.assertRaisesRegex(ValueError,'Mixed identified/anonymous'):self.prep(self.fixture().replace('ID=a;',''))
    def test_collision_refused(self):
        with self.assertRaisesRegex(ValueError,'collision'):
            self.prep(self.fixture()+'s\tTiberius\texon\t1\t9\t.\t+\t.\tID=t.CDS;Parent=t\n')
    def test_mismatched_strands_refused(self):
        with self.assertRaisesRegex(ValueError,'sequence/strand mismatch'):
            self.prep(self.fixture().replace('CDS\t20\t28\t.\t+','CDS\t20\t28\t.\t-'))
    def test_shared_cds_split_before_grouping(self):
        text=self.fixture().replace('Parent=t','Parent=t,u')
        text+='s\tTiberius\tmRNA\t1\t100\t.\t+\t.\tID=u;Parent=g\n'
        out,r=self.prep(text)
        self.assertEqual(out.count('ID=t.CDS;'),2);self.assertEqual(out.count('ID=u.CDS;'),2)
        self.assertEqual(r['shared_features_split'],2);self.assertNotIn('Parent=t,u',out)
    def test_shared_internal_feature_refused(self):
        text=self.fixture().replace('ID=t;Parent=g','ID=t;Parent=g,h')+'s\tTiberius\tgene\t1\t100\t.\t+\t.\tID=h\n'
        with self.assertRaisesRegex(ValueError,'internal or referenced'):self.prep(text)
    def test_external_transcript_identifier_preserved(self):
        text=self.fixture().replace('ID=t;Parent=g','ID=t;Parent=g;transcript_id=RNA123.1').replace('ID=a;Parent=t','ID=a;Parent=t;transcript_id=RNA123.1').replace('ID=b;Parent=t','ID=b;Parent=t;transcript_id=RNA123.1')
        out,r=self.prep(text);self.assertIn('original_transcript_id=RNA123.1',out)
    def test_syntax_override_rejects_mismatch(self):
        with self.assertRaisesRegex(ValueError,'input-format gtf'):self.prep(input_format='gtf')
    def test_mixed_syntax_refused(self):
        with self.assertRaisesRegex(ValueError,'mixed GTF/GFF3'):
            self.prep(self.fixture()+'s\tUnknown\texon\t1\t9\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
    def test_phase_normalized_for_unknown_source(self):
        out,r=self.prep(self.fixture('Unknown')+'s\tUnknown\texon\t1\t9\t.\t+\t2\tID=e;Parent=t\n')
        self.assertIn('original_phase=2',out);self.assertIn('\texon\t1\t9\t.\t+\t.\t',out)
    def test_phase_conflict_on_excluded_organelle_does_not_run(self):
        metadata=self.root/'ncbi.json';metadata.write_text(json.dumps({'reports':[record('s'),record('mt','Mitochondrion')]}))
        text=self.fixture('Unknown')+MT+'mt\ttest\texon\t1\t90\t.\t+\t1\tID=em;Parent=tm;original_phase=2\n'
        out,r=self.prep(text,remove_organelles=True,sequence_metadata=str(metadata))
        self.assertNotIn('\nmt\t',out);self.assertEqual(r['organelle_filter']['removed_feature_rows'],4)
    def test_provenance_independent_of_normalization_mode(self):
        self.source.write_text('##gff-version 3\n'+self.fixture('Tiberius'))
        settings=dict(input=str(self.source),assembly_accession=ACC,locus_tag_prefix='TEST')
        a=g.resolve_settings({**settings,'normalization':'standard'});b=g.resolve_settings({**settings,'normalization':'strict'})
        self.assertEqual(a['inference_by_source'],b['inference_by_source'])
    def test_invalid_policy_rejected(self):
        with self.assertRaisesRegex(ValueError,'grouping overrides'):self.prep(normalization='strict',cds_grouping='preserve')
    def test_direct_gene_cds_existing_group_preserved(self):
        text='s\tUnknown\tgene\t1\t100\t.\t+\t.\tID=g\ns\tUnknown\tCDS\t1\t90\t.\t+\t0\tID=protein;Parent=g\n'
        out,r=self.prep(text);self.assertIn('ID=protein;Parent=g',out);self.assertNotIn('\tmRNA\t',out)
    def test_unknown_gtf_without_tool_hint(self):
        text='s\tUnknown\tgene\t1\t100\t.\t+\t.\tgene_id "g";\n'
        text+='s\tUnknown\ttranscript\t1\t100\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'
        text+='s\tUnknown\tCDS\t1\t90\t.\t+\t0\tgene_id "g"; transcript_id "t";\n'
        out,r=self.prep(text,input_format='gtf');self.assertIn('ID=t;Parent=g',out);self.assertIn('\tmRNA\t',out)
    def test_shared_anonymous_leaves_expanded(self):
        text=self.fixture()+'s\tTiberius\tmRNA\t1\t100\t.\t+\t.\tID=u;Parent=g\ns\tTiberius\texon\t1\t9\t.\t+\t.\tParent=t,u\n'
        out,r=self.prep(text);self.assertNotIn('Parent=t,u',out);self.assertEqual(r['shared_features_split'],1)
    def test_strict_still_produces_single_parent_rows(self):
        text=self.fixture().replace('\ttranscript\t','\tmRNA\t')+'s\tTiberius\tmRNA\t1\t100\t.\t+\t.\tID=u;Parent=g\ns\tTiberius\texon\t1\t9\t.\t+\t.\tID=e;Parent=t,u\n'
        out,r=self.prep(text,normalization='strict');self.assertNotIn('Parent=t,u',out);self.assertIn('ID=e.ena_parent1',out)
    def test_shared_output_idempotent(self):
        text=self.fixture().replace('Parent=t','Parent=t,u')+'s\tTiberius\tmRNA\t1\t100\t.\t+\t.\tID=u;Parent=g\n'
        out,r=self.prep(text)
        # prep adds declarations, so feed only features on the second pass.
        again,r=self.prep('\n'.join(line for line in out.splitlines() if not line.startswith('#'))+'\n')
        self.assertEqual(out,again)
