import uuid
from unittest import TestCase

from tgmi.transcripts import TranscriptDBWriter

from refseqdb import refseq
import os


class RefSeqReader(TestCase):

    def setUp(self):

        ncbi_mappings_NM_000059 = [{
            'chrom': '13',
            'strand': '+',
            'exonStarts': [
                32889616, 32890558, 32893213, 32899212, 32900237, 32900378, 32900635, 32903579,
                32905055, 32906408, 32910401, 32918694, 32920963, 32928997, 32930564, 32931878,
                32936659, 32937315, 32944538, 32945092, 32950806, 32953453, 32953886, 32954143,
                32968825, 32971034, 32972298
            ],
            'exonEnds': [
                32889804, 32890664, 32893462, 32899321, 32900287, 32900419, 32900750, 32903629,
                32905167, 32907524, 32915333, 32918790, 32921033, 32929425, 32930746, 32932066,
                32936830, 32937670, 32944694, 32945237, 32950928, 32953652, 32954050, 32954282,
                32969070, 32971181, 32973809
            ],
            'coding_start': 32890597,
            'coding_end': 32972906,
            'exon_cigars': [
                '188=', '106=', '249=', '109=', '50=', '41=', '115=', '50=', '112=', '1116=', '4932=', '96=', '70=',
                '389=1X38=', '182=', '188=', '171=', '355=', '156=', '145=', '122=', '199=', '164=', '139=', '245=',
                '147=', '1511='
            ]
        }]

        self.mappings = { 'NM_000059': ncbi_mappings_NM_000059}

        self.outfn = str(uuid.uuid4())
        self.out_incl = open(self.outfn + '_included.txt', 'w')
        self.out_excl = open(self.outfn + '_excluded.txt', 'w')

        # Initialize tdb_writer
        columns = [
            'ID',
            'VERSION',
            'HGNC_ID',
            'INFO',
            'STRAND',
            'CHROM',
            'START',
            'END',
            'EXONS',
            'CODING_START',
            'CODING_END',
            'SEQUENCE',
            'CDNA_CODING_START',
            'CDNA_CODING_END'
        ]
        self.tdb_writer = TranscriptDBWriter(self.outfn, columns=columns)


    def tearDown(self):
        os.remove(self.outfn + '_included.txt')
        os.remove(self.outfn + '_excluded.txt')


    def test_process_record_included(self):

        record = []
        for line in open('test/unit/data/NM_000059_interim.txt'):
            line = line.strip()
            record.append(line)

        refseq.process_record(self.mappings, self.tdb_writer, self.out_incl, self.out_excl, record)

        r = self.tdb_writer._records['13'][0]

        assert r[:6] == ['NM_000059', '3', 'HGNC:1101', '+/84.2kb/27/11.4kb/3418', '+', '13']
        assert r[6] == 32889616
        assert r[7] == 32973809

        assert r[8].split(',')[2] == '32893213-32893462'
        assert r[8].split(',')[-1] == '32972298-32973809'

        assert r[9] == 32890597
        assert r[10] == 32972906

        assert r[11][:10] == 'gtggcgcgag'.upper()
        assert r[11][-16:] == 'caaattggcactgatt'.upper()

        assert r[12] == 228
        assert r[13] == 10484

        self.out_incl.close()
        self.out_excl.close()

        for line in open(self.outfn + '_included.txt'):
            line = line.strip()
            break
        assert line == 'NM_000059\tHGNC:1101'


    def test_process_features_section(self):

        features_section = []
        for line in open('test/unit/data/features.txt'):
            line = line.strip()
            features_section.append(line)

        cdna_coding_start, cdna_coding_end, hgncid = refseq.process_features_section(features_section)
        assert cdna_coding_start == '228'
        assert cdna_coding_end == '10484'
        assert hgncid == 'HGNC:1101'

        features_section[10] = 'CDS             join(168..260,262..741)'
        cdna_coding_start, cdna_coding_end, hgncid = refseq.process_features_section(features_section)
        assert cdna_coding_start is None

        features_section = ['FEATURES']
        cdna_coding_start, cdna_coding_end, hgncid = refseq.process_features_section(features_section)
        assert cdna_coding_start == ''
        assert cdna_coding_end == ''
        assert hgncid == ''


    def test_process_origin_section(self):
        pass


    def test_load_data_into_transcript_object(self):
        pass

    def test_check_for_issues(self):
        pass

    def test_split_sections(self):
        record = []
        for line in open('test/unit/data/NM_000059_interim.txt'):
            line = line.strip()
            record.append(line)

        result = refseq.split_sections(record)
        assert len(result['LOCUS']) == 1
        assert len(result['VERSION']) == 1
        assert len(result['FEATURES']) == 290
        assert len(result['ORIGIN']) == 192