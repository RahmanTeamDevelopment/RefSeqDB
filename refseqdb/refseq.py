from __future__ import division
from tgmi.transcripts import Transcript, Exon
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


def process_refseq_file(fn, mappings, tdb_writer, out_incl, out_excl):
    """Process RefSeq data file"""
    record = []
    for line in gzip.open(fn, 'r'):
        line = line.strip()
        if line.startswith('LOCUS'):
            process_record(mappings, tdb_writer, out_incl, out_excl, record)
            record = []
        record.append(line)
    process_record(mappings, tdb_writer, out_incl, out_excl, record)


def process_record(mappings, tdb_writer, out_incl, out_excl, record):
    """Process a RefSeq record"""

    if len(record) == 0:
        return

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

            if line.startswith('CDS') and '..' in line:
                coords = line.split()[1]
                if coords.startswith('join'):
                    out_excl.write('\t'.join([id, 'joined_cds']) + '\n')
                    return
                tmp = coords.split('..')
                cdna_coding_start, cdna_coding_end = tmp[0], tmp[-1]

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
    if id not in mappings:
        out_excl.write('\t'.join([transcript.id, 'no_mapping']) + '\n')
        return

    # Multiple UCSC mapping
    if len(mappings[id]) > 1:
        out_excl.write('\t'.join([transcript.id, 'multiple_mapping']) + '\n')
        return

    # Single mapping
    mapping = mappings[id][0]

    # Set chrom, strand, exons, start and end of transcript
    transcript.chrom = mapping['chrom']
    transcript.strand = mapping['strand']

    # Set comma-delimited exon cigar strings
    transcript.exon_cigars = ','.join(mapping['exon_cigars']) if mapping['exon_cigars'] is not None else '.'

    transcript.exons = []
    for i in range(len(mapping['exonStarts'])):
        transcript.exons.append(Exon(str(mapping['exonStarts'][i]) + '-' + str(mapping['exonEnds'][i])))

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
