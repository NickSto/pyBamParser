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

    def test_basic(self):
        coord_pairs = ((0, None), (1, 781), (2, 782), (284, 1064), (285, None))
        self._test_many_coords(781, '284M', 163, 284, coord_pairs)

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
        coord_pairs = ((1, 6307), (2, 6306), (3, 6305), (286, 6022))
        self._test_many_coords(6022, '286M', 83, 286, coord_pairs)

    def test_reverse_insertion(self):
        coord_pairs = ((1, 8299), (25, 8275), (26, None), (34, None), (35, 8274))
        self._test_many_coords(8027, '248M9I25M', 83, 282, coord_pairs)

    def test_reverse_deletion(self):
        coord_pairs = ((1, 3123), (2, 3122), (17, 3107), (18, 3105))
        self._test_many_coords(2840, '266M1D17M', 83, 283, coord_pairs)

    def test_reverse_deletion_toy(self):
        coord_pairs = ((1, 1300), (100, 1201), (101, 1100), (200, 1001))
        self._test_many_coords(1001, '100M100D100M', 83, 200, coord_pairs)

    # The following are taken from Figure 1 of Li et al. 2009 which introduced the SAM format.
    def test_Li_r001p(self):
        coord_pairs = ((1, 7), (8, 14), (9, None), (10, None), (11, 15), (14, 18), (15, 20))
        self._test_many_coords(7, '8M2I4M1D3M', 163, 17, coord_pairs)

    def test_Li_r002p(self):
        coord_pairs = ((1, None), (3, None), (4, 9), (9, 14), (10, None), (11, 15), (14, 18))
        self._test_many_coords(9, '3S6M1P1I4M', 0, 14, coord_pairs)

    def test_Li_r003p(self):
        coord_pairs = ((0, None), (1, 9), (6, 14), (7, None))
        self._test_many_coords(9, '5H6M', 0, 6, coord_pairs)

    def test_Li_r004p(self):
        coord_pairs = ((1, 16), (6, 21), (7, 36), (11, 40))
        self._test_many_coords(16, '6M14N5M', 0, 11, coord_pairs)

    def test_Li_r003n(self):
        coord_pairs = ((1, 33), (5, 29))
        self._test_many_coords(29, '6H5M', 16, 5, coord_pairs)

    def test_Li_r001n(self):
        coord_pairs = ((1, 45), (9, 37))
        self._test_many_coords(37, '9M', 83, 9, coord_pairs)


def logging_setup():
    logging.basicConfig(stream=sys.stderr, level=logging.WARNING, format='%(message)s')
    for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
        level_name = logging.getLevelName(level)
        logging.addLevelName(level, level_name.capitalize())

if __name__ == '__main__':
  main(sys.argv)
