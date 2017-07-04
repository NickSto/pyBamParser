#!/usr/bin/env python2
import sys
import logging
import argparse
import unittest
import pyBamParser.read
import cigar

ARG_DEFAULTS = {'log':sys.stderr, 'volume':logging.WARNING}
DESCRIPTION = """"""


def make_argparser():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.set_defaults(**ARG_DEFAULTS)
    parser.add_argument('-l', '--log', type=argparse.FileType('w'),
        help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
    parser.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL)
    parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
    parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)
    return parser


def main(argv):
    parser = make_argparser()
    args = parser.parse_args(argv[1:])
    logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
    for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
        level_name = logging.getLevelName(level)
        logging.addLevelName(level, level_name.capitalize())
    unittest.main()

"""
# Find an example of a type of CIGAR string:
$ samtools view duplex.down10.bam | awk '! and($2, 16) && $3 == "chrM" && $6 ~ /^[0-9]*M[0-9]*S$/ {print $6}' | sort | uniq -c | sort -g -k 1 | head
# Get the rest of the info for the alignment:
$ samtools view duplex.down10.bam | awkt '! and($2, 16) && $3 == "chrM" && $6 == "255M12S"' | cut -f -8
# Check the name is unique:
$ samtools view duplex.down10.bam | grep AACCCAAGTGACTGCTCGCACTTA | wc -l
# Try converting some coordinates:
$ samtools view duplex.down10.bam | grep AACCCAAGTGACTGCTCGCACTTA | ./cigar.py 255
"""

class CigarConversionTest(unittest.TestCase):

    def _test_many_coords(self, pos, cigar_str, flags, readlen, coord_pairs):
        read = pyBamParser.read.BAMRead('12345678901234567890123456789012', None)
        cigar_list = cigar.split_cigar(cigar_str)
        blocks = read._get_contiguous_blocks(pos, cigar_list, flags & 16, readlen)
        for read_coord, ref_coord in coord_pairs:
            result = read._to_ref_coord(blocks, read_coord)
            try:
                self.assertEqual(ref_coord, result)
            except AssertionError:
                logging.warn('Failed {}: {} -> {} (got {} instead)'
                             .format(cigar_str, read_coord, ref_coord, result))
                raise

    def test_insertion(self):
        coord_pairs = ((159, 8270), (160, None), (168, None), (169, 8271))
        self._test_many_coords(8112, '159M9I115M', 99, 283, coord_pairs)

    def test_deletion(self):
        coord_pairs = ((111, 3105), (112, 3107))
        self._test_many_coords(2995, '111M1D172M', 99, 283, coord_pairs)

    def test_left_soft_padding(self):
        coord_pairs = ((11, None), (12, 5059), (13, 5060))
        self._test_many_coords(5059, '11S267M', 99, 278, coord_pairs)

    def test_right_soft_padding(self):
        coord_pairs = ((255, 6528), (256, None))
        self._test_many_coords(6274, '255M12S', 163, 267, coord_pairs)

    def test_reverse(self):
        coord_pairs = ((1, 6307), (2, 6306), (3, 6305), (284, 6022))
        self._test_many_coords(6022, '286M', 83, 286, coord_pairs)

    def test_reverse_insertion(self):
        coord_pairs = ((1, 8299), (25, 8275), (26, None), (34, None), (35, 8274))
        self._test_many_coords(8027, '248M9I25M', 83, 282, coord_pairs)

    def test_reverse_deletion(self):
        coord_pairs = ((1, 3123), (2, 3122), (17, 3107), (18, 3105))
        self._test_many_coords(2840, '266M1D17M', 83, 283, coord_pairs)


def logging_setup():
    logging.basicConfig(stream=sys.stderr, level=logging.WARNING, format='%(message)s')
    for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
        level_name = logging.getLevelName(level)
        logging.addLevelName(level, level_name.capitalize())

if __name__ == '__main__':
  main(sys.argv)
