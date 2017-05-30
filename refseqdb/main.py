from __future__ import division
from utils.transcripts import TranscriptDBWriter, Transcript, Exon
from urllib2 import urlopen
import gzip
import os
import sys
from ftplib import FTP

########################################################################################################################

# Return list of available RefSeq data files
def access_refseq_files():
    ret = []
    ftp = FTP('ftp.ncbi.nlm.nih.gov', timeout = 3600)
    ftp.login()
    ftp.cwd('/refseq/H_sapiens/mRNA_Prot/')
    files = []
    for x in ftp.nlst():
        if x.startswith('human.') and x.endswith('.rna.gbff.gz'):
            files.append(x)
    i = 0
    while True:
        i += 1
        fn = 'human.'+str(i)+'.rna.gbff.gz'
        if fn in files: ret.append('ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/'+fn)
        else: break
    return ret

# Download and read UCSC mappings data
def read_ucsc_mapping(build):
    ret = dict()

    # Download appropriate data file depending on build
    if build == 'GRCh38': url = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz'
    else: url = 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz'
    f = urlopen(url)
    with open('ucsc.gz', "wb") as datafile: datafile.write(f.read())

    # Iterate through the lines of data file
    for line in gzip.open('ucsc.gz'):
        line = line.strip()
        if line == '': continue
        cols = line.split()
        if not cols[1].startswith('NM_'): continue

        chrom = cols[2][3:]
        allowed = [str(x) for x in range(1,23)]+['X','Y','MT']
        if chrom not in allowed: continue

        key = cols[1]
        if key not in ret.keys(): ret[key] = []

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

        ret[key].append(mapping)

    os.remove('ucsc.gz')
    return ret

# Process a single RefSeq data file
def process_refseq_file(UCSC_mappings, tdb_writer, out_incl, out_excl):
    record = []
    for line in gzip.open('refseqdata.gz', 'r'):
        line = line.strip()
        if line.startswith('LOCUS'):
            process_record(UCSC_mappings, tdb_writer, out_incl, out_excl, record)
            record = []
        record.append(line)
    process_record(UCSC_mappings, tdb_writer, out_incl, out_excl, record)

# Process a single record
def process_record(UCSC_mappings, tdb_writer, out_incl, out_excl, record):
    if len(record) == 0: return

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
    ID = locusline.split()[1]

    if not ID.startswith('NM_'): return

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
        if line.startswith('ORIGIN'): continue
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

    # Output nothing in case of missing data
    if cdna_coding_start == '' or cdna_coding_end == '' or sequence == '' or version == '' or hgncid == '': return

    # Initialize transcript and set ID, version and hgnc_id
    transcript = Transcript()
    transcript.id = ID
    transcript.version = version
    transcript.hgnc_id = hgncid

    transcript.sequence = sequence
    transcript.cdna_coding_start = cdna_coding_start
    transcript.cdna_coding_end = cdna_coding_end

    # No UCSC mapping
    if ID not in UCSC_mappings:
        out_excl.write('\t'.join([transcript.id, 'no_mapping']) + '\n')
        return

    # Multiple UCSC mapping
    if len(UCSC_mappings[ID]) > 1:
        out_excl.write('\t'.join([transcript.id, 'multiple_mapping']) + '\n')
        return

    # Single mapping
    mapping = UCSC_mappings[ID][0]

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

    # Add transcript to database
    tdb_writer.add(transcript)
    out_incl.write('\t'.join([transcript.id, transcript.hgnc_id])+'\n')


# Split record to LOCUS, VERSION, FEATURES and ORIGIN sections
def split_sections(record):
    ret = {'LOCUS': [], 'VERSION': [], 'FEATURES': [], 'ORIGIN': []}
    keywords = ['LOCUS', 'DEFINITION', 'ACCESSION', 'VERSION', 'KEYWORDS', 'SOURCE', 'REFERENCE', 'FEATURES', 'ORIGIN']
    key = ''
    for line in record:
        for k in keywords:
            if line.startswith(k):
                key = k
                break
        if key in ret.keys(): ret[key].append(line)
    return ret

# Initialize progress information
def init_progress_info():
    sys.stdout.write('\rProcessing transcripts ... 0.0%')
    sys.stdout.flush()

# Print out progress information
def print_progress_info(counter, N):
    x = round(100 * counter / N, 1)
    x = min(x, 100.0)
    sys.stdout.write('\rProcessing transcripts ... ' + str(x) + '%')
    sys.stdout.flush()

# Finalize progress information
def finalize_progress_info():
    sys.stdout.write('\rProcessing transcripts ... 100.0%')
    sys.stdout.flush()
    print ' - Done.'


def main(options):

    ver = 'x.x.x'

    columns = ['ID', 'VERSION', 'HGNC_ID', 'INFO', 'STRAND', 'CHROM', 'START', 'END', 'EXONS', 'CODING_START', 'CODING_END',
               'SEQUENCE', 'CDNA_CODING_START', 'CDNA_CODING_END']

    tdb_writer = TranscriptDBWriter(options.output, source='refseq_db ' + ver, build=options.build, columns=columns)

    out_incl = open(options.output + '_included.txt', 'w')
    out_excl = open(options.output + '_excluded.txt', 'w')

    # Retrieve list of available RefSeq data files
    sys.stdout.write('\rAccessing RefSeq data files ... ')
    sys.stdout.flush()
    urls = access_refseq_files()
    print '- Done.'

    # Download and read UCSC mapping data
    sys.stdout.write('\rDownloading and reading UCSC mapping data ... ')
    sys.stdout.flush()
    UCSC_mappings = read_ucsc_mapping(options.build)
    print '- Done.'
    print ''

    # Initialize progress info
    init_progress_info()

    unknown_errors = []
    # Iterate through available RefSeq data files
    for i in range(len(urls)):

        # Download RefSeq data file and avoid timeout error
        while True:
            try:
                f = urlopen(urls[i])
                break
            except:
                pass
        with open('refseqdata.gz', "wb") as datafile:
            datafile.write(f.read())

        # Process RefSeq data file
        print_progress_info(i + 1, len(urls))
        process_refseq_file(UCSC_mappings, tdb_writer, out_incl, out_excl)

        # Remove RefSeq data file
        os.remove('refseqdata.gz')

    # Finalize progress info
    finalize_progress_info()

    # Finalize transcript database
    tdb_writer.finalize()

    out_incl.close()

    # Goodbye message
    print ''
    if len(unknown_errors) > 0: print 'Unknown errors: ' + str(unknown_errors) + '\n'


