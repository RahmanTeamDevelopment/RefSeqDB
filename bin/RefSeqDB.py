#!env/bin/python

import datetime
from optparse import OptionParser
from refseqdb.main import main
import sys

# Version
ver = '0.2.1'

# Command line argument parsing
descr = 'RefSeqDB v'+ver
parser = OptionParser(usage='RefSeqDB/env/bin/refseqdb <options>', version=ver, description=descr)
#parser.add_option('-b', default='GRCh37', dest='build', action='store', help="Genome build [default value: %default]")
parser.add_option('-m', default='ncbi', dest='mapping', action='store', help="Mapping to use (ncbi or ucsc) [default value: %default]")
parser.add_option('-o', default='output', dest='output', action='store', help="Output file name prefix [default value: %default]")
(options, args) = parser.parse_args()

if options.mapping not in ['ucsc', 'ncbi']:
    print '\nAllowed values for option -m are \"ucsc\" or \"ncbi\".\n'
    sys.exit()

# Welcome message
print '\n' + '='*100
now = str(datetime.datetime.now())
now = now[:now.find('.')]
print 'RefSeqDB v{} started: {}\n'.format(ver, now)

if options.mapping == 'ncbi':
    print '- RefSeq interim release, 2017-01-13 (incl. mapping) -\n'
else:
    print '- Latest RefSeq release + UCSC mapping -\n'

main(ver, options)

print '\nOutput files created: '
print ' - {}.gz (+ .tbi)'.format(options.output)
print ' - {}_included.txt'.format(options.output)
print ' - {}_excluded.txt'.format(options.output)

now = str(datetime.datetime.now())
now = now[:now.find('.')]
print '\nFinished: {}'.format(now)
print '='*100 + '\n'
