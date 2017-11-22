from urllib2 import urlopen
import gzip
import pysam


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
            'exonEnds': [x for x in cols[10].split(',') if x != ''],
            'exon_cigars': None
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


def download_ncbi_mapping(build):
    """Download NCBI mappings data"""

    if build=='GRCh37':
        ftpsite = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/GRCh37.p13_interim_annotation/'
        files = [
            'interim_GRCh37.p13_knownrefseq_alignments_2017-01-13.bam',
            'interim_GRCh37.p13_top_level_2017-01-13.gff3.gz'
        ]
    else:
        ftpsite = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/GRCh38.p10_interim_annotation/'
        files = [
            'interim_GRCh38.p10_knownrefseq_alignments_2017-01-13.bam',
            'interim_GRCh38.p10_top_level_2017-01-13.gff3.gz'
        ]

    for fn in files:
        f = urlopen(ftpsite + fn)
        with open(fn, "wb") as datafile:
            datafile.write(f.read())


def read_ncbi_mapping(bam_fn, gff3_fn):
    """Read NCBI mappings data"""

    ret = {}

    # Read strand, coding_start and coding_end from GFF3 file
    gff3_data = process_gff3_file(gff3_fn)

    for line in pysam.view(bam_fn):
        line = line.strip()
        cols = line.split()

        id = cols[0]
        if '.' in id:
            id = id[:id.find('.')]

        if not id.startswith('NM_'):
            continue

        # Extract chromosome name
        chrom_nc = cols[2]
        if not chrom_nc.startswith('NC_'):
            continue
        if '.' in chrom_nc:
            chrom_nc = chrom_nc[:chrom_nc.find('.')]
        chrom = translate_chromosome_name_from_nc_id(chrom_nc)

        start_pos = int(cols[3]) - 1

        # Process cigar string and calculate exon coordinates
        cigar = cols[5]
        cigar_list = split_cigar(cigar)
        exon_lengths, intron_lengths, exon_cigars = break_into_exons(cigar_list)
        exon_starts, exon_ends = calculate_exon_coordinates(exon_lengths, intron_lengths, start_pos)

        # Extract strand, coding_start, coding_end from the GFF3 data
        # Note: the followings are not correct in case of multiple mappings, however transcripts with multiple mappings
        # are excluded anyway. If transcripts with multiple mappings are included in the future, the code should match the
        # corresponding mappings read from the two different sources and not always add data from the 0th mapping.
        strand = gff3_data[id][0][0]
        coding_start = gff3_data[id][0][1]
        coding_end = gff3_data[id][0][2]

        # Reverse order of exon coordinates for reverse-stranded transcripts
        if strand == '-':
            exon_starts = exon_starts[::-1]
            exon_ends = exon_ends[::-1]
            exon_cigars = exon_cigars[::-1]

        mapping = {
            'chrom': chrom,
            'strand': strand,
            'exonStarts': exon_starts,
            'exonEnds': exon_ends,
            'coding_start': coding_start,
            'coding_end': coding_end,
            'exon_cigars': exon_cigars
        }

        if id not in ret:
            ret[id] = []
        ret[id].append(mapping)

    return ret


def process_gff3_file(gff3_fn):
    """Read strand, coding_start and coding_end from GFF3 file"""

    ret = {}

    reading_id = ''

    reading = False
    for line in gzip.open(gff3_fn):
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

        if cols[2] not in ['mRNA', 'CDS']:
            continue

        if cols[2] == 'mRNA':

            if reading:
                if id not in ret:
                    ret[id] = []
                coding_start, coding_end = coding_endpoints(strand, min_cds, max_cds)
                ret[id].append((strand, coding_start, coding_end))
                reading = False

            id = extract_id_from_info_field(cols[8])

            if id.startswith('NM_'):
                min_cds = None
                max_cds = None
                strand = cols[6]
                reading = True
                reading_id = extract_gff3_record_id_from_info_field(cols[8])

            continue

        parent = extract_parent_from_info_field(cols[8])
        if parent != reading_id:
            continue

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
        coding_start, coding_end = coding_endpoints(strand, min_cds, max_cds)
        ret[id].append((strand, coding_start, coding_end))

    return ret


def coding_endpoints(strand, min_cds, max_cds):
    """Determine coding_start and coding_end coordinates"""

    if strand == '+':
        return min_cds - 1, max_cds - 1
    else:
        return max_cds - 1, min_cds - 1


def extract_id_from_info_field(info):
    """Helper function to extract transcript ID from info string"""

    for tag in info.split(';'):
        if tag.startswith('Dbxref='):
            for x in tag.split(','):
                if x.startswith('Genbank:'):
                    id = x[8:]
                    if '.' in id:
                        id = id[:id.find('.')]
                    return id


def extract_parent_from_info_field(info):
    """Helper function to extract parent ID from info string"""

    for tag in info.split(';'):
        if tag.startswith('Parent='):
            return tag[7:]


def extract_gff3_record_id_from_info_field(info):
    """Helper function to extract GFF3 record ID from info string"""

    for tag in info.split(';'):
        if tag.startswith('ID='):
            return tag[3:]


def split_cigar(cig):
    """Split cigar string to list"""

    ret = []
    s = ''
    for c in cig:
        s += c
        if c in ['N', '=', 'X', 'I', 'D']:
            ret.append(s)
            s = ''
        cig = cig[1:]
    return ret


def break_into_exons(cigar_list):
    """Break cigar list to exons"""

    exons = []
    intron_lengths = []
    e = []
    for x in cigar_list:
        if x[-1] == 'N':
            exons.append(e)
            intron_lengths.append(int(x[:-1]))
            e = []
        else:
            e.append(x)
    exons.append(e)

    exon_lengths = []
    for e in exons:
        total = 0
        for x in e:
            if x[-1] in ['=', 'X', 'D']:
                total += int(x[:-1])
        exon_lengths.append(total)

    exon_cigars = []
    for e in exons:
        exon_cigar = ''
        for x in e:
            exon_cigar += x
        exon_cigars.append(exon_cigar)

    return exon_lengths, intron_lengths, exon_cigars


def calculate_exon_coordinates(exon_lengths, intron_lengths, start_pos):
    """Calculate exon coordinates"""

    exon_starts = []
    exon_ends = []

    p = start_pos
    for i in range(len(exon_lengths)):
        exon_starts.append(p)
        exon_ends.append(p + exon_lengths[i])
        if i == len(exon_lengths) - 1:
            break
        p += exon_lengths[i] + intron_lengths[i]

    return exon_starts, exon_ends


def translate_chromosome_name_from_nc_id(nc):
    """Helper function to translate chromosome name from NC ID"""

    if nc == 'NC_000023':
        return 'X'
    if nc == 'NC_000024':
        return 'Y'
    if nc == 'NC_012920':
        return 'MT'
    return str(int(nc[3:]))


