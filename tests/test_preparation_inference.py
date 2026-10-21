import unittest
import test_gffprep as base
m=base.m
class PreparationInferenceTests(unittest.TestCase):
    setUp=base.Tests.setUp
    runprep=base.Tests.runprep
    def test_inference(self):
        _,report,out=self.runprep(inference='ab initio prediction:ANNEVO:2.2.3')
        self.assertEqual(out.count('inference=ab initio prediction:ANNEVO:2.2.3'),2)
        self.assertEqual(report['submission_metadata']['CDS_inference'],'ab initio prediction:ANNEVO:2.2.3')
        _,_,again=self.runprep(out,inference='ab initio prediction:ANNEVO:2.2.3')
        self.assertEqual(out,again)
        with self.assertRaisesRegex(m.Invalid,'existing inference differs'):
            self.runprep(out,inference='ab initio prediction:Other:1')
    def test_bad_inference(self):
        for value in ['', ' ', 'ab\ninitio', 'ab\té']:
            with self.assertRaisesRegex(m.Invalid,'inference'):self.runprep(inference=value)
