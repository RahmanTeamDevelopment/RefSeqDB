from urllib2 import urlopen
import gzip


def download_ucsc_mapping(build, fn):
    """Download UCSC mappings data"""

    # Download appropriate data file depending on build
    if build == 'GRCh38':
        url = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz'
    else:
        url = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz'
    f = urlopen(url)
    with open(fn, "wb") as datafile:
        datafile.write(f.read())


def download_ncbi_mapping(fn):
    """Download NCBI mappings data"""

    f = urlopen('ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/GRCh37.p13_interim_annotation/interim_GRCh37.p13_top_level_2017-01-13.gff3.gz')
    with open(fn, "wb") as datafile:
        datafile.write(f.read())


def read_ucsc_mapping(fn):
    """Read UCSC mappings data"""

    ret = dict()

    allowed_chroms = [str(x) for x in range(1, 23)] + ['X', 'Y', 'MT']

    # Iterate through the lines of data file
    for line in gzip.open(fn):
        line = line.strip()
        if line == '':
            continue
        cols = line.split()

        if not cols[1].startswith('NM_'):
            continue

        chrom = cols[2][3:]

        if chrom not in allowed_chroms:
            continue

        key = cols[1]

        # Save mapping in a dict
        mapping = {
            'chrom': chrom,
            'strand': cols[3],
            'coding_start': cols[6],
            'coding_end': cols[7],
            'exonStarts': [x for x in cols[9].split(',') if x != ''],
            'exonEnds': [x for x in cols[10].split(',') if x != '']
        }
        if mapping['strand'] == '-':
            mapping['exonStarts'] = mapping['exonStarts'][::-1]
            mapping['exonEnds'] = mapping['exonEnds'][::-1]
        if mapping['strand'] == '+':
            mapping['coding_start'] = int(cols[6])
            mapping['coding_end'] = int(cols[7]) - 1
        else:
            mapping['coding_start'] = int(cols[7]) - 1
            mapping['coding_end'] = int(cols[6])

        # Add possible mapping to list
        try:
            ret[key].append(mapping)
        except:
            ret[key] = [mapping]

    return ret


def read_ncbi_mapping(fn):
    """Read NCBI mappings data"""

    ret = {}
    exon_starts = []
    exon_ends = []
    reading_id = ''

    reading = False
    for line in gzip.open(fn):
        line = line.strip()
        if line == '' or line.startswith('#'):
            continue

        cols = line.split('\t')

        if not cols[0].startswith('NC_'):
            continue

        '''
        if cols[1] != 'BestRefSeq':
            continue
        '''

        if cols[2] not in ['mRNA', 'exon', 'CDS']:
            continue

        if cols[2] == 'mRNA':

            if reading:
                if id not in ret:
                    ret[id] = []
                ret[id].append(finalize_mapping(strand, min_cds, max_cds, chrom, exon_starts, exon_ends))
                reading = False

            id = extract_id(cols[8])

            if id.startswith('NM_'):
                exon_starts = []
                exon_ends = []
                min_cds = None
                max_cds = None
                strand = cols[6]
                chrom_nc = cols[0]
                chrom_nc = chrom_nc[:chrom_nc.find('.')]
                chrom = translate_chromosome_name_from_nc_id(chrom_nc)
                reading = True
                reading_id = extract_gff3_record_id(cols[8])

            continue

        parent = extract_parent(cols[8])
        if parent != reading_id:
            continue

        if reading and cols[2] == 'exon':
            exon_starts.append(int(cols[3])-1)
            exon_ends.append(int(cols[4]))

        if reading and cols[2] == 'CDS':
            if min_cds is None:
                min_cds = int(cols[3])
                max_cds = int(cols[4])
            else:
                if int(cols[3]) < min_cds:
                    min_cds = int(cols[3])
                if int(cols[4]) > max_cds:
                    max_cds = int(cols[4])

    if reading:
        if id not in ret:
            ret[id] = []
        ret[id].append(finalize_mapping(strand, min_cds, max_cds, chrom, exon_starts, exon_ends))

    return ret


def extract_id(info):
    """Helper function to extract transcript ID from info string"""

    for tag in info.split(';'):
        if tag.startswith('Dbxref='):
            for x in tag.split(','):
                if x.startswith('Genbank:'):
                    id = x[8:]
                    if '.' in id:
                        id = id[:id.find('.')]
                    return id


def extract_parent(info):
    """Helper function to extract parent ID from info string"""

    for tag in info.split(';'):
        if tag.startswith('Parent='):
            return tag[7:]


def extract_gff3_record_id(info):
    """Helper function to extract GFF3 record ID from info string"""

    for tag in info.split(';'):
        if tag.startswith('ID='):
            return tag[3:]


def translate_chromosome_name_from_nc_id(nc):
    """Helper function to translate chromosome name from NC ID"""

    if nc == 'NC_000023':
        return 'X'
    if nc == 'NC_000024':
        return 'Y'
    if nc == 'NC_012920':
        return 'MT'
    return str(int(nc[3:]))


def finalize_mapping(strand, min_cds, max_cds, chrom, exon_starts, exon_ends):
    """Helper function to finalize mapping data"""

    if strand == '+':
        coding_start = min_cds - 1
        coding_end = max_cds - 1
    else:
        coding_start = max_cds - 1
        coding_end = min_cds - 1

    return {
            'chrom': chrom,
            'strand': strand,
            'exonStarts': exon_starts,
            'exonEnds': exon_ends,
            'coding_start': coding_start,
            'coding_end': coding_end
    }