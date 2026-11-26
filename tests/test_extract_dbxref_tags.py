import unittest
from tools.extract_dbxref_tags import extract


def row(label, example):
    return f'<tr><td>{label}</td><td>Description</td><td>{example}</td><td>&nbsp;</td></tr>'


def table(*rows):
    return '<table>' + ''.join(rows) + '</table>'


class ExtractTests(unittest.TestCase):
    def test_markup_and_case(self):
        tags, issues = extract(table(row('<a>Demo(in host)</a>&nbsp;', '/db_xref=&ldquo;Demo(in host):1&rdquo;')))
        self.assertEqual(tags, {'Demo(in host)'})
        self.assertEqual(issues, [])

    def test_disagreement(self):
        tags, issues = extract(table(row('Demo', 'db_xref=DEMO:1<br>db_xref=Demo:2')))
        self.assertEqual(tags, {'DEMO', 'Demo'})
        self.assertEqual(issues, [('Demo', ['DEMO', 'Demo'])])

    def test_duplicate(self):
        r = row('Demo', 'db_xref=Demo:1')
        with self.assertRaises(ValueError):
            extract(table(r, r))

    def test_missing_example(self):
        with self.assertRaises(ValueError):
            extract(table(row('Demo', 'db_xref=Demo:1'), row('Other', 'missing')))

    def test_ambiguous_table(self):
        t = table(row('Demo', 'db_xref=Demo:1'))
        with self.assertRaises(ValueError):
            extract(t + t)

    def test_empty(self):
        with self.assertRaises(ValueError):
            extract('<table></table>')

    def test_truncated(self):
        with self.assertRaises(ValueError):
            extract(table(row('Demo', 'db_xref=Demo:1'))[:-8])
