
from unittest import TestCase
from refseqdb import mapping


class TestNCBIMapping(TestCase):

    def test_read_ncbi_mapping(self):

        self.ncbi_mappings_NM_000059 = [{
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

        self.ncbi_mappings_NM_001039348 = [{
            'chrom': '2',
            'strand': '-',
            'exonStarts': [
                56150845, 56150033, 56149494, 56145353, 56144799, 56108746, 56104880, 56103757,
                56102080, 56098134, 56097854, 56093096
            ],
            'exonEnds': [
                56151298, 56150074, 56149582, 56145402, 56145186, 56108869, 56105000, 56103877,
                56102200, 56098258, 56098050, 56094369
            ],
            'coding_start': 56149574,
            'coding_end': 56094207,
            'exon_cigars': ['453=', '41=', '88=', '49=', '387=', '123=', '120=', '120=', '120=', '124=', '196=',
                            '1273=']
        }]

        self.ncbi_mappings_forward = mapping.read_ncbi_mapping(
            'test/unit/data/interim_test_forward.bam',
            'test/unit/data/interim_test_forward.gff3.gz'
        )

        self.ncbi_mappings_reverse = mapping.read_ncbi_mapping(
            'test/unit/data/interim_test_reverse.bam',
            'test/unit/data/interim_test_reverse.gff3.gz'
        )

        assert set(self.ncbi_mappings_forward.keys()) == {'NM_000059'}
        assert set(self.ncbi_mappings_reverse.keys()) == {'NM_001039348','NM_001039349'}

        assert self.ncbi_mappings_forward['NM_000059'] == self.ncbi_mappings_NM_000059
        assert self.ncbi_mappings_reverse['NM_001039348'] == self.ncbi_mappings_NM_001039348

    def test_process_gff3_file(self):

        result = mapping.process_gff3_file('test/unit/data/interim_test_forward.gff3.gz')
        assert set(result.keys()) == {'NM_000059'}
        assert result['NM_000059'][0] == ('+', 32890597, 32972906)

        result = mapping.process_gff3_file('test/unit/data/interim_test_reverse.gff3.gz')
        assert set(result.keys()) == {'NM_001039348','NM_001039349'}
        assert result['NM_001039348'][0] == ('-', 56149574, 56094207)

    def test_extract_id_from_info_field(self):
        info = 'ID=rna7031;Parent=gene4539;Dbxref=GeneID:2202,Genbank:NM_001039348.2,HGNC:HGNC:3218,MIM:601548;Name=NM_001039348.2;gbkey=mRNA;gene=EFEMP1;product=EGF containing fibulin like extracellular matrix protein 1%2C transcript variant 2;transcript_id=NM_001039348.2'
        assert mapping.extract_id_from_info_field(info) == 'NM_001039348'

    def test_extract_parent_from_info_field(self):
        info = 'ID=cds5242;Parent=rna7031;Dbxref=CCDS:CCDS1857.1,GeneID:2202,Genbank:NP_001034437.1,HGNC:HGNC:3218,MIM:601548;Name=NP_001034437.1;gbkey=CDS;gene=EFEMP1;product=EGF-containing fibulin-like extracellular matrix protein 1 precursor;protein_id=NP_001034437.1'
        assert mapping.extract_parent_from_info_field(info) == 'rna7031'

    def test_extract_gff3_record_id_from_info_field(self):
        info = 'ID=rna7031;Parent=gene4539;Dbxref=GeneID:2202,Genbank:NM_001039348.2,HGNC:HGNC:3218,MIM:601548;Name=NM_001039348.2;gbkey=mRNA;gene=EFEMP1;product=EGF containing fibulin like extracellular matrix protein 1%2C transcript variant 2;transcript_id=NM_001039348.2'
        assert mapping.extract_gff3_record_id_from_info_field(info) == 'rna7031'

    def test_split_cigar(self):
        assert mapping.split_cigar('10=20N15=1D8=12I6X') == ['10=', '20N','15=', '1D', '8=', '12I', '6X']

    def test_break_into_exons(self):

        cig = mapping.split_cigar('188=754N106=2549N249=')
        assert mapping.break_into_exons(cig)[0] == [188, 106, 249]
        assert mapping.break_into_exons(cig)[1] == [754, 2549]
        assert mapping.break_into_exons(cig)[2] == ['188=', '106=', '249=']

        cig = mapping.split_cigar('5X188=754N106=2549N249=3X10=')
        assert mapping.break_into_exons(cig)[0] == [193, 106, 262]
        assert mapping.break_into_exons(cig)[1] == [754, 2549]
        assert mapping.break_into_exons(cig)[2] == ['5X188=', '106=', '249=3X10=']

        cig = mapping.split_cigar('188=754N106=7I2549N249=')
        assert mapping.break_into_exons(cig)[0] == [188, 106, 249]
        assert mapping.break_into_exons(cig)[1] == [754, 2549]
        assert mapping.break_into_exons(cig)[2] == ['188=', '106=7I', '249=']

        cig = mapping.split_cigar('188=754N106=7D2549N249=')
        assert mapping.break_into_exons(cig)[0] == [188, 113, 249]
        assert mapping.break_into_exons(cig)[1] == [754, 2549]
        assert mapping.break_into_exons(cig)[2] == ['188=', '106=7D', '249=']

    def test_calculate_exon_coordinates(self):
        result = mapping.calculate_exon_coordinates([3, 4], [6], 10)
        assert result[0] == [10, 19]
        assert result[1] == [13, 23]

        result = mapping.calculate_exon_coordinates([20, 30, 15], [12, 7], 1000)
        assert result[0] == [1000, 1032, 1069]
        assert result[1] == [1020, 1062, 1084]

    def test_translate_chromosome_name_from_nc_id(self):
        assert mapping.translate_chromosome_name_from_nc_id('NC_000004') == '4'
        assert mapping.translate_chromosome_name_from_nc_id('NC_000019') == '19'
        assert mapping.translate_chromosome_name_from_nc_id('NC_000023') == 'X'
        assert mapping.translate_chromosome_name_from_nc_id('NC_000024') == 'Y'
        assert mapping.translate_chromosome_name_from_nc_id('NC_012920') == 'MT'

