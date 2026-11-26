import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
import test_gffprep as base
import ena_gffprep as cli

class OptionalTaxonomyTests(unittest.TestCase):
    setUp=base.Tests.setUp
    runprep=base.Tests.runprep

    def test_no_taxon_no_species(self):
        _,r,out=self.runprep(taxon_id=None)
        self.assertNotIn('##species',out)
        self.assertIsNone(r['submission_metadata']['taxon_id'])
        self.assertEqual(r['submission_metadata']['taxonomy_source'],'not_provided')

    def test_preserve_input_species_when_no_argument(self):
        directive='##species https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/9606'
        _,r,out=self.runprep(directive+'\n'+base.BASE,taxon_id=None)
        self.assertEqual(out.count('##species'),1)
        self.assertIn(directive,out)
        self.assertEqual(r['submission_metadata']['taxon_id'],9606)
        self.assertEqual(r['submission_metadata']['taxonomy_source'],'input_header')

    def test_preserve_input_species_with_annevo_adapter(self):
        text=base.BASE.replace('ID=c1','ID=t1.CDS.1',1)
        text=text.replace('ID=c1','ID=t1.CDS.2')
        _,_,out=self.runprep('##species https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/\n'+text,
                             source_format='annevo',taxon_id=None)
        self.assertIn('##species https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/',out)

    def test_explicit_taxon_still_adds_directive(self):
        _,r,out=self.runprep(taxon_id=9606)
        self.assertIn('##species https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/',out)
        self.assertEqual(r['submission_metadata']['taxonomy_source'],'explicit')

    def test_invalid_taxon_and_input_conflicts(self):
        for value in [0,-1,True]:
            with self.assertRaisesRegex(base.m.Invalid,'positive integer'):self.runprep(taxon_id=value)
        text='##species https://www.ncbi.nlm.nih.gov/datasets/taxonomy/9606/\n'+base.BASE
        with self.assertRaisesRegex(base.m.Invalid,'conflicts'):self.runprep(text,taxon_id=10090)
        text='##species https://www.ncbi.nlm.nih.gov/datasets/taxonomy/10090/\n'+text
        with self.assertRaisesRegex(base.m.Invalid,'Conflicting input'):self.runprep(text,taxon_id=None)

    def test_batch_without_taxon_column(self):
        (self.d/'a.gff').write_text(base.BASE)
        (self.d/'batch.tsv').write_text('input\tassembly_accession\tlocus_tag_prefix\n'
            'a.gff\tGCA_123456789.1\tTEST\n')
        with contextlib.redirect_stdout(io.StringIO()):
            code=cli.batch_prepare(self.d/'batch.tsv',{'bioproject':'PRJEB123456','no_inference':True,'out_root':str(self.d/'out')},self.d/'summary.json')
        self.assertEqual(code,0)
        path=next((self.d/'out').glob('*/preparation.json'))
        self.assertIsNone(json.loads(path.read_text())['submission_metadata']['taxon_id'])

    def test_cli_without_taxon(self):
        (self.d/'a.gff').write_text(base.BASE)
        with contextlib.redirect_stdout(io.StringIO()):
            code=cli.main(['--input',str(self.d/'a.gff'),'--assembly-accession','GCA_123456789.1','--bioproject','PRJEB123456',
                '--locus-tag-prefix','TEST','--no-inference','--output-dir',str(self.d/'out')])
        self.assertEqual(code,0)
        self.assertTrue((self.d/'out/manifest.txt').is_file())
