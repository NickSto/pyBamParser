#!/usr/bin/env python
import os
import sys
from pyBamParser.read import BAMRead
from pyBamParser.bam import Reader
from optparse import OptionParser

def main():

  FUNCTIONS = {
    'BAMRead.get_indels':BAMRead_get_indels,
    'BAMRead.indel_at':BAMRead_indel_at,
  }

  OPT_DEFAULTS = {'str':'', 'int':0, 'bool':False}
  USAGE = "USAGE: %prog [options] function.to.test reads.bam"
  DESCRIPTION = """Run test on a given function and input BAM and print results.
  Give one of the following function names: """+', '.join(FUNCTIONS)
  EPILOG = """ """

  parser = OptionParser(usage=USAGE, description=DESCRIPTION, epilog=EPILOG)

  # parser.add_option('-s', '--str', dest='str',
  #   default=OPT_DEFAULTS.get('str'), help='default: %default')
  # parser.add_option('-i', '--int', dest='int', type='int',
  #   default=OPT_DEFAULTS.get('int'), help='')
  # parser.add_option('-b', '--bool', dest='bool', action='store_const',
  #   const=not OPT_DEFAULTS.get('bool'), default=OPT_DEFAULTS.get('bool'),
  #   help='')

  (options, arguments) = parser.parse_args()

  if len(arguments) == 2:
    (function, bamfilename) = arguments
  else:
    parser.print_help()
    fail('Provide a function name and a BAM file.')

  if function not in FUNCTIONS:
    fail('Error: function "'+function+'" not supported. Please pick one from '
      +'the list: '+', '.join(FUNCTIONS))
  if not os.path.exists(bamfilename):
    fail('Error: cannot find BAM file "'+bamfilename+'"')

  bam_reader = Reader(bamfilename)

  FUNCTIONS[function](bam_reader)


def BAMRead_get_indels(bam_reader):
  for read in bam_reader:
    print "\t".join([read.get_read_name(),str(read.get_position()),read.get_sam_cigar()])
    # print "contains indel at "+str(2280)+": "+str(read.indel_at(2280))
    (insertions, deletions) = read.get_indels()
    if insertions:
      print "\t"+str(insertions)
    if deletions:
      print "\t"+str(deletions)


def BAMRead_indel_at(bam_reader):
  pass


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()