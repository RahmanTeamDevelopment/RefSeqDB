#!env/bin/python

import datetime
from optparse import OptionParser
from refseqdb.main import main

# Version
ver = 'x.x'

# Command line argument parsing
descr = 'refseq_db v'+ver
parser = OptionParser(usage='python path/to/refseq_db/refseq_db.py <options>', version=ver, description=descr)
parser.add_option('-b', "--build", default='GRCh37', dest='build', action='store', help="Genome build [default value: %default]")
parser.add_option('-o', "--out", default='output.txt', dest='output', action='store', help="Output file name [default value: %default]")
(options, args) = parser.parse_args()

# Welcome message
print '\n'+'='*100
print 'refseq_db v'+ver+' started: '+str(datetime.datetime.now())+'\n'

main(options)


print 'Finished: '+str(datetime.datetime.now())
print '='*100+'\n'
