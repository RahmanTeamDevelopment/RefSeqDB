from __future__ import division
from tgmi.transcripts import TranscriptDBWriter
from urllib2 import urlopen
import os
import sys
import refseq
import mapping


def init_progress_info():
    """Initialize progress information"""

    sys.stdout.write('\rProcessing transcripts ... 0.0%')
    sys.stdout.flush()


def print_progress_info(counter, N):
    """Print out progress information"""

    x = round(100 * counter / N, 1)
    x = min(x, 100.0)
    sys.stdout.write('\rProcessing transcripts ... ' + str(x) + '%')
    sys.stdout.flush()


def finalize_progress_info():
    """Finalize progress information"""

    sys.stdout.write('\rProcessing transcripts ... 100.0%')
    sys.stdout.flush()
    print ' - Done.'


def main(ver, options):
    """Main function"""

    columns = ['ID', 'VERSION', 'HGNC_ID', 'INFO', 'STRAND', 'CHROM', 'START', 'END', 'EXONS', 'CODING_START',
               'CODING_END', 'SEQUENCE', 'CDNA_CODING_START', 'CDNA_CODING_END']

    tdb_writer = TranscriptDBWriter(options.output, source='refseq_db ' + ver, build='GRCh37', columns=columns)

    # Initialize output files and write headers
    out_incl = open(options.output + '_included.txt', 'w')
    out_incl.write('\t'.join(['#ID', 'HGNCID']) + '\n')
    out_excl = open(options.output + '_excluded.txt', 'w')
    out_excl.write('\t'.join(['#ID', 'Reason']) + '\n')

    # Retrieve list of available RefSeq data files
    sys.stdout.write('\rAccessing RefSeq data files ... ')
    sys.stdout.flush()

    if options.mapping == 'ucsc':
        urls = refseq.access_refseq_files()
    else:
        urls = ['ftp://ftp.ncbi.nlm.nih.gov/genomes/Homo_sapiens/GRCh37.p13_interim_annotation/interim_GRCh37.p13_rna.gbk.gz']
    print '- Done.'

    # Download and read UCSC mapping data
    sys.stdout.write('\rDownloading and reading mapping data ... ')
    sys.stdout.flush()
    if options.mapping == 'ucsc':
        mapping.download_ucsc_mapping('GRCh37', 'ucsc_mapping.gz')
        mappings = mapping.read_ucsc_mapping('ucsc_mapping.gz')
        os.remove('ucsc_mapping.gz')
    else:
        mapping.download_ncbi_mapping('ncbi_mapping.gz')
        mappings = mapping.read_ncbi_mapping('ncbi_mapping.gz')
        os.remove('ncbi_mapping.gz')
    print '- Done.'
    print ''

    # Initialize progress info
    init_progress_info()

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
        refseq.process_refseq_file('refseqdata.gz', mappings, tdb_writer, out_incl, out_excl)

        # Remove RefSeq data file
        os.remove('refseqdata.gz')

    # Finalize progress info
    finalize_progress_info()

    # Finalize transcript database
    tdb_writer.finalize()

    # Close output files
    out_incl.close()
    out_excl.close()
