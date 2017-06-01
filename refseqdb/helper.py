from __future__ import division
from utils.transcripts import Transcript, Exon
from urllib2 import urlopen
from ftplib import FTP
import gzip


def access_refseq_files():
    """Return list of available RefSeq data files"""

    ret = []
    ftp = FTP('ftp.ncbi.nlm.nih.gov', timeout=3600)
    ftp.login()
    ftp.cwd('/refseq/H_sapiens/mRNA_Prot/')
    files = []
    for x in ftp.nlst():
        if x.startswith('human.') and x.endswith('.rna.gbff.gz'):
            files.append(x)
    i = 0
    while True:
        i += 1
        fn = 'human.' + str(i) + '.rna.gbff.gz'
        if fn in files:
            ret.append('ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/' + fn)
        else:
            break
    return ret


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


def process_refseq_file(ucsc_mappings, tdb_writer, out_incl, out_excl):
    """Process RefSeq data file"""
    record = []
    for line in gzip.open('refseqdata.gz', 'r'):
        line = line.strip()
        if line.startswith('LOCUS'):
            process_record(ucsc_mappings, tdb_writer, out_incl, out_excl, record)
            record = []
        record.append(line)
    process_record(ucsc_mappings, tdb_writer, out_incl, out_excl, record)


def process_record(ucsc_mappings, tdb_writer, out_incl, out_excl, record):
    """Process a RefSeq record"""

    if len(record) == 0:
        return

    cdna_exons = []
    dna = ''
    cdna_coding_start = ''
    cdna_coding_end = ''
    hgncid = ''
    chromosome = ''

    # Split record to LOCUS, VERSION, FEATURES and ORIGIN sections
    sections = split_sections(record)

    # Retrieve transcript ID (NM)
    locusline = sections['LOCUS'][0]
    id = locusline.split()[1]

    if not id.startswith('NM_'):
        return

    # Retrieve version
    versionline = sections['VERSION'][0]
    version = versionline.split()[1]
    version = version.split('.')[1]

    # Scan FEATURES section to retrieve chromosome, exons, CDS and HGNCID
    for line in sections['FEATURES']:
        try:
            if line.startswith('/chromosome=') and chromosome == '':
                chromosome = line[13:-1]

            if line.startswith('exon'):
                coords = line.split()[1]
                [start, end] = coords.split('..')
                cdna_exons.append(start + '-' + end)

            if line.startswith('CDS'):
                coords = line.split()[1]
                cdna_coding_start, cdna_coding_end = coords.split('..')

            if line.startswith('/db_xref="HGNC:') and hgncid == '':
                hgncid = line[15:-1]
        except:
            pass

    # Retrieve transcript sequence
    for line in sections['ORIGIN']:
        if line.startswith('ORIGIN'):
            continue
        dna += line
    sequence = ''
    for x in dna:
        if x in ['a', 'c', 'g', 't']:
            sequence += x.upper()
    sequence_trimmed = ''
    for cdna_exon in cdna_exons:
        [start, end] = cdna_exon.split('-')
        sequence_trimmed += sequence[int(start) - 1:int(end)]
    sequence = sequence_trimmed

    # Missing CDS
    if cdna_coding_start == '' or cdna_coding_end == '':
        out_excl.write('\t'.join([id, 'missing_cds']) + '\n')
        return

    # Missing sequence
    if sequence == '':
        out_excl.write('\t'.join([id, 'missing_sequence']) + '\n')
        return

    # Missing version
    if version == '':
        out_excl.write('\t'.join([id, 'missing_version']) + '\n')
        return

    # Missing HGNC ID
    if hgncid == '':
        out_excl.write('\t'.join([id, 'missing_hgncid']) + '\n')
        return

    # Initialize transcript and set ID, version and hgnc_id
    transcript = Transcript(id=id, version=version, hgnc_id=hgncid)
    transcript.sequence = sequence
    transcript.cdna_coding_start = cdna_coding_start
    transcript.cdna_coding_end = cdna_coding_end

    # No UCSC mapping
    if id not in ucsc_mappings:
        out_excl.write('\t'.join([transcript.id, 'no_mapping']) + '\n')
        return

    # Multiple UCSC mapping
    if len(ucsc_mappings[id]) > 1:
        out_excl.write('\t'.join([transcript.id, 'multiple_mapping']) + '\n')
        return

    # Single mapping
    mapping = ucsc_mappings[id][0]

    # Set chrom, strand, exons, start and end of transcript
    transcript.chrom = mapping['chrom']
    transcript.strand = mapping['strand']

    transcript.exons = []
    for i in range(len(mapping['exonStarts'])):
        transcript.exons.append(Exon(mapping['exonStarts'][i] + '-' + mapping['exonEnds'][i]))

    if transcript.strand == '+':
        transcript.start = transcript.exons[0].start
        transcript.end = transcript.exons[-1].end
    else:
        transcript.start = transcript.exons[-1].start
        transcript.end = transcript.exons[0].end

    transcript.coding_start = int(mapping['coding_start'])
    transcript.coding_end = int(mapping['coding_end'])

    # Finalize transcript
    transcript.finalize()

    # Incorrect CDS length
    if transcript.get_cds_length() % 3 != 0:
        out_excl.write('\t'.join([id, 'incorrect_cds_length']) + '\n')
        return

    # Add transcript to database
    tdb_writer.add(transcript)
    out_incl.write('\t'.join([transcript.id, transcript.hgnc_id]) + '\n')


def split_sections(record):
    """Split record to LOCUS, VERSION, FEATURES and ORIGIN sections"""

    ret = {'LOCUS': [], 'VERSION': [], 'FEATURES': [], 'ORIGIN': []}
    keywords = ['LOCUS', 'DEFINITION', 'ACCESSION', 'VERSION', 'KEYWORDS', 'SOURCE', 'REFERENCE', 'FEATURES', 'ORIGIN']
    key = ''
    for line in record:
        for k in keywords:
            if line.startswith(k):
                key = k
                break
        if key in ret.keys():
            ret[key].append(line)
    return ret

