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
    read_name      = read.get_read_name()
    flag           = read.get_flag()
    reference      = read.get_reference()
    reference_name = read.get_reference_name()
    reference_id   = read.get_reference_id()
    rnext          = read.get_rnext()
    rnext_name     = read.get_rnext_name()
    position       = read.get_position()
    end_position   = read.get_end_position()
    pnext          = read.get_pnext()
    mapq           = read.get_mapq()
    sam_cigar      = read.get_sam_cigar()
    cigar          = read.get_cigar()
    t_len          = read.get_t_len()
    seq            = read.get_seq()
    l_seq          = read.get_l_seq()
    qual           = read.get_qual()
    sam_qual       = read.get_sam_qual()
    qual_tuple     = read.get_qual_tuple()
    is_seq_revcomp = read.is_seq_reverse_complement()
    read_group     = read.get_read_group()
    sam_aux        = read.get_sam_aux()
    print "%-14s: %-14s '%s'" % ("read_name", type(read_name), read_name)
    print "%-14s: %-14s '%s'" % ("flag", type(flag), flag)
    print "%-14s: %-14s '%s'" % ("reference", type(reference), reference)
    print "%-14s: %-14s '%s'" % ("reference_name", type(reference_name), reference_name)
    print "%-14s: %-14s '%s'" % ("reference_id", type(reference_id), reference_id)
    print "%-14s: %-14s '%s'" % ("rnext", type(rnext), rnext)
    print "%-14s: %-14s '%s'" % ("rnext_name", type(rnext_name), rnext_name)
    print "%-14s: %-14s '%s'" % ("position", type(position), position)
    print "%-14s: %-14s '%s'" % ("end_position", type(end_position), end_position)
    print "%-14s: %-14s '%s'" % ("pnext", type(pnext), pnext)
    print "%-14s: %-14s '%s'" % ("mapq", type(mapq), mapq)
    print "%-14s: %-14s '%s'" % ("sam_cigar", type(sam_cigar), sam_cigar)
    print "%-14s: %-14s '%s'" % ("cigar", type(cigar), cigar)
    print "%-14s: %-14s '%s'" % ("t_len", type(t_len), t_len)
    print "%-14s: %-14s '%s'" % ("seq", type(seq), seq)
    print "%-14s: %-14s '%s'" % ("l_seq", type(l_seq), l_seq)
    print "%-14s: %-14s '%s'" % ("qual", type(qual), qual)
    print "%-14s: %-14s '%s'" % ("sam_qual", type(sam_qual), sam_qual)
    print "%-14s: %-14s '%s'" % ("qual_tuple", type(qual_tuple), qual_tuple)
    print "%-14s: %-14s '%s'" % ("is_seq_revcomp", type(is_seq_revcomp), is_seq_revcomp)
    print "%-14s: %-14s '%s'" % ("read_group", type(read_group), read_group)
    print "%-14s: %-14s '%s'" % ("sam_aux", type(sam_aux), sam_aux)
    print
    # print "contains indel at "+str(2280)+": "+str(read.indel_at(2280))
    # (insertions, deletions) = read.get_indels()
    # if insertions:
    #   print "\t"+str(insertions)
    # if deletions:
    #   print "\t"+str(deletions)


def fail(message):
  sys.stderr.write(message+"\n")
  sys.exit(1)

if __name__ == "__main__":
  main()