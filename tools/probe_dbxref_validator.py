"""Generate synthetic standalone-validator probes; no external data required."""
import argparse
from pathlib import Path

CASES = {
    '00_baseline': None,
    '01_geneid': 'GeneID:12345',
    '02_bold_upper': 'BOLD:TEST123',
    '03_bold_mixed': 'Bold:TEST123',
    '04_interpro_upper': 'InterPro:IPR000001',
    '05_interpro_lower': 'Interpro:IPR000001',
    '06_phytozome': 'Phytozome:TEST123',
    '07_phytozme': 'Phytozme:TEST123',
    '08_unknown_database': 'NonexistentDatabase:12345',
    '09_invalid_syntax': 'InvalidReference',
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    out = parser.parse_args().output
    out.mkdir(parents=True, exist_ok=False)
    for name, value in CASES.items():
        attrs = 'ID=gene1;locus_tag=TEST_00001'
        if value is not None:
            attrs += ';Dbxref=' + value
        lines = [
            '##gff-version 3', '#!insdc-gff-version 0.5.1',
            '##sequence-region TESTSEQ 1 1000',
            f'TESTSEQ\ttest\tgene\t1\t300\t.\t+\t.\t{attrs}',
            'TESTSEQ\ttest\tmRNA\t1\t300\t.\t+\t.\tID=mrna1;Parent=gene1',
            'TESTSEQ\ttest\texon\t1\t300\t.\t+\t.\tParent=mrna1',
            'TESTSEQ\ttest\tCDS\t1\t300\t.\t+\t0\tID=cds1;Parent=mrna1;product=hypothetical protein;transl_table=1',
        ]
        (out / (name + '.gff3')).write_text('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
