import tempfile,unittest,sys,gzip
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import ena_add_inference as m
class InferenceTests(unittest.TestCase):
 def test_provenance_only_and_failure_cleanup(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);source=p/'a.gff';dest=p/'b.gff.gz'
   text='##gff-version 3\ns\tx\tCDS\t1\t60\t.\t+\t0\tID=c;Parent=t\n'
   source.write_text(text);self.assertEqual(m.add(source,dest,'ab initio prediction:ANNEVO:2.2.3'),1)
   with gzip.open(dest,'rt') as f:result=f.read()
   self.assertEqual(result.replace(';inference=ab initio prediction:ANNEVO:2.2.3',''),text)
   self.assertEqual(source.read_text(),text)
   with self.assertRaisesRegex(ValueError,'already exists'):m.add(source,dest,'x')
   source.write_text(text.replace('ID=c','inference=old;ID=c'))
   with self.assertRaisesRegex(ValueError,'Existing inference'):m.add(source,p/'fail.gz','x')
   self.assertFalse((p/'fail.gz').exists())
