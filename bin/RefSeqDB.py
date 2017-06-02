#!env/bin/python

import datetime
from optparse import OptionParser
from refseqdb.main import main

# Version
ver = '0.1.0'

# Command line argument parsing
descr = 'RefSeqDB v'+ver
parser = OptionParser(version=ver, description=descr)
parser.add_option('-b', "--build", default='GRCh37', dest='build', action='store', help="Genome build [default value: %default]")
parser.add_option('-o', "--out", default='output', dest='output', action='store', help="Output file name prefix [default value: %default]")
(options, args) = parser.parse_args()

# Welcome message
print '\n'+'='*100
print 'RefSeqDB v'+ver+' started: '+str(datetime.datetime.now())+'\n'

main(ver, options)

print '\nFinished: '+str(datetime.datetime.now())
print '='*100+'\n'
