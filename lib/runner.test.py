#!/usr/bin/env python
import os
import sys
from pyBamParser.read import BAMRead
from pyBamParser.bam import Reader
from optparse import OptionParser

OPT_DEFAULTS = {'str':'', 'int':0, 'float':0.0, 'bool':False}
USAGE = "USAGE: %prog [options]"
DESCRIPTION = """"""
EPILOG = """"""

def main():

  parser = OptionParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)

  parser.add_option('-s', '--str', dest='str',
    default=OPT_DEFAULTS.get('str'), help='default: %default')
  parser.add_option('-i', '--int', dest='int', type='int',
    default=OPT_DEFAULTS.get('int'), help='')
  parser.add_option('-f', '--float', dest='float', type='float',
    default=OPT_DEFAULTS.get('float'), help='')
  parser.add_option('-b', '--bool', dest='bool', action='store_const',
    const=not OPT_DEFAULTS.get('bool'), default=OPT_DEFAULTS.get('bool'),
    help='')

  (options, arguments) = parser.parse_args()

  if not arguments:
    parser.print_help()
    fail('Give a BAM file')
  else:
    bamfilename = arguments[0]

  bam_reader = Reader(bamfilename)
  for read in bam_reader:
    print "\t".join([read.get_read_name(),str(read.get_position()),read.get_sam_cigar()])
    # print "contains indel at "+str(2280)+": "+str(read.indel_at(2280))
    (insertions, deletions) = read.get_indels()
    if insertions:
      print "\t"+str(insertions)
    if deletions:
      print "\t"+str(deletions)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()