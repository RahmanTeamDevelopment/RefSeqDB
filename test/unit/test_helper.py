"""Unit tests for helper.py"""


from unittest import TestCase
from utils.transcripts import TranscriptDBWriter
from refseqdb import helper
import uuid
import os


class TestHelper(TestCase):

    def setUp(self):

        # Initialize output files
        self.outfn = str(uuid.uuid4())
        self.out_incl = open(self.outfn + '_included.txt', 'w')
        self.out_excl = open(self.outfn + '_excluded.txt', 'w')

        # UCSC mapping of NM_000059:
        self.ucsc_mappings_NM_000059 = {
            'NM_000059': [{
                'chrom': '13',
                'strand': '+',
                'exonStarts': [
                    '32889616', '32890558', '32893213', '32899212', '32900237', '32900378', '32900635', '32903579',
                    '32905055', '32906408', '32910401', '32918694', '32920963', '32928997', '32930564', '32931878',
                    '32936659', '32937315', '32944538', '32945092', '32950806', '32953453', '32953886', '32954143',
                    '32968825', '32971034', '32972298'
                ],
                'exonEnds': [
                    '32889804', '32890664', '32893462', '32899321', '32900287', '32900419', '32900750', '32903629',
                    '32905167', '32907524', '32915333', '32918790', '32921033', '32929425', '32930746', '32932066',
                    '32936830', '32937670', '32944694', '32945237', '32950928', '32953652', '32954050', '32954282',
                    '32969070', '32971181', '32973809'],
                'coding_start': 32890597,
                'coding_end': 32972906
            }]
        }

        self.ucsc_mappings = helper.read_ucsc_mapping('test/unit/data/ucsc_mapping.gz')

        # Initialize tdb_writer
        columns = ['ID', 'VERSION', 'HGNC_ID', 'INFO', 'STRAND', 'CHROM', 'START', 'END', 'EXONS', 'CODING_START',
                   'CODING_END', 'SEQUENCE', 'CDNA_CODING_START', 'CDNA_CODING_END']
        self.tdb_writer = TranscriptDBWriter(self.outfn, columns=columns)


    def tearDown(self):
        # Remove output files
        os.remove(self.outfn + '_included.txt')
        os.remove(self.outfn + '_excluded.txt')


    def test_split_sections(self):
        record = []
        for line in open('test/unit/data/NM_000059.txt'):
            line = line.strip()
            record.append(line)

        result = helper.split_sections(record)
        assert len(result['LOCUS']) == 1
        assert len(result['VERSION']) == 1
        assert len(result['FEATURES']) == 478
        assert len(result['ORIGIN']) == 192


    def test_read_ucsc_mapping(self):
        assert sum(map(len,[self.ucsc_mappings[k] for k in self.ucsc_mappings])) == 48216
        assert self.ucsc_mappings['NM_000059'][0] == self.ucsc_mappings_NM_000059['NM_000059'][0]


    def test_process_record_included(self):

        # NM_000059 (included):

        record = []
        for line in open('test/unit/data/NM_000059.txt'):
            line = line.strip()
            record.append(line)

        helper.process_record(self.ucsc_mappings_NM_000059, self.tdb_writer, self.out_incl, self.out_excl, record)

        r = self.tdb_writer._records['13'][0]

        assert r[:6] == ['NM_000059', '3', 'HGNC:1101', '+/84.2kb/27/11.4kb/3418', '+', '13']
        assert r[6] == 32889616
        assert r[7] == 32973809
        assert r[9] == 32890597
        assert r[10] == 32972906

        assert r[8].split(',')[2] == '32893213-32893462'
        assert r[8].split(',')[-1] == '32972298-32973809'

        self.out_incl.close()
        self.out_excl.close()

        for line in open(self.outfn+'_included.txt'):
            line = line.strip()
            break
        assert line == 'NM_000059\tHGNC:1101'


    def test_process_record_incomplete(self):

        # NM_001301851 (missing hgncid):

        record = []
        for line in open('test/unit/data/NM_001301851.txt'):
            line = line.strip()
            record.append(line)

        helper.process_record(self.ucsc_mappings, self.tdb_writer, self.out_incl, self.out_excl, record)

        self.out_incl.close()
        self.out_excl.close()

        for line in open(self.outfn + '_excluded.txt'):
            line = line.strip()
            break
        assert line == 'NM_001301851\tmissing_hgncid'


    def test_process_record_multiple_mapping(self):

        # NM_152585 (multiple mapping)

        record = []
        for line in open('test/unit/data/NM_152585.txt'):
            line = line.strip()
            record.append(line)

        helper.process_record(self.ucsc_mappings, self.tdb_writer, self.out_incl, self.out_excl, record)

        self.out_incl.close()
        self.out_excl.close()

        for line in open(self.outfn + '_excluded.txt'):
            line = line.strip()
            break
        assert line == 'NM_152585\tmultiple_mapping'


