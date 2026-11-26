import unittest
from ena_gffprep import parser, resolve_settings

class BioProjectTests(unittest.TestCase):
    def test_cli_flag(self):
        self.assertEqual(parser().parse_args(['--bioproject','PRJEB123456']).bioproject, 'PRJEB123456')

    def test_setting_recognized_in_config_and_batch(self):
        from ena_gffprep import SETTING_KEYS
        self.assertIn('bioproject', SETTING_KEYS)

if __name__ == '__main__':
    unittest.main()
