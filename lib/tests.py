#!/usr/bin/env python
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
        logging.info(blocks)
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


class CigarGetIndelsTest(unittest.TestCase):

    def _test_read(self, pos, cigar_str, flags, readlen, expected):
        read = pyBamParser.read.BAMRead('12345678901234567890123456789012', None)
        cigar_list = cigar.split_cigar(cigar_str)
        blocks = read._get_contiguous_blocks(pos, cigar_list, flags & 16, readlen)
        logging.info(blocks)
        indels = read._get_indels(blocks, flags & 16)
        self.assertEqual(indels, expected)

    def test_M(self):
        indels = ([], [])
        self._test_read(1, '251M', 99, 251, indels)

    def test_SM(self):
        indels = ([], [])
        self._test_read(1, '3S248M', 99, 251, indels)

    def test_MS(self):
        indels = ([], [])
        self._test_read(111, '205M46S', 163, 251, indels)

    def test_MIMa(self):
        indels = ([(5, 2)], [])
        self._test_read(2, '4M2I245M', 147, 251, indels)

    def test_MIMb(self):
        indels = ([(310, 1)], [])
        self._test_read(199, '112M1I138M', 83, 251, indels)

    def test_MDMa(self):
        indels = ([], [(3106, 1)])
        self._test_read(2941, '166M1D85M', 99, 251, indels)

    def test_MDMb(self):
        indels = ([], [(1900, 1)])
        self._test_read(1785, '116M1D135M', 83, 251, indels)

    def test_SMS(self):
        indels = ([], [])
        self._test_read(1, '10S156M85S', 163, 251, indels)

    def test_SMDM(self):
        indels = ([], [(2528, 1)])
        self._test_read(2526, '11S3M1D237M', 83, 251, indels)

    def test_SMIM(self):
        indels = ([(556, 3)], [])
        self._test_read(554, '38S3M3I207M', 83, 251, indels)

    def test_MDMS(self):
        indels = ([], [(3106, 1)])
        self._test_read(2883, '224M1D26M1S', 163, 251, indels)

    def test_MIMS(self):
        indels = ([(4269, 1)], [])
        self._test_read(4099, '171M1I11M68S', 99, 251, indels)

    def test_MIMDM(self):
        indels = ([(14880, 2)], [(14883, 2)])
        self._test_read(14640, '241M2I3M2D5M', 163, 251, indels)

    def test_MIMIM(self):
        indels = ([(16573, 2), (16569, 39)], [])
        self._test_read(16365, '205M39I4M2I1M', 83, 251, indels)

    def test_MIMDMSa(self):
        indels = ([(10172, 1)], [(10175, 2)])
        self._test_read(9931, '242M1I3M2D2M3S', 99, 251, indels)

    def test_MIMDMSb(self):
        indels = ([(16568, 62)], [(16570, 2)])
        self._test_read(16390, '179M62I2M2D2M6S', 99, 251, indels)

    def test_SMDMDM(self):
        indels = ([], [(6808, 1), (6806, 1)])
        self._test_read(6800, '3S7M1D1M1D240M', 147, 251, indels)

    def test_SMIMIM(self):
        indels = ([(5980, 1), (5979, 1)], [])
        self._test_read(5975, '10S5M1I1M1I233M', 83, 251, indels)

    def test_MIMIMS(self):
        indels = ([(11348, 3), (11354, 2)], [])
        self._test_read(11127, '222M3I6M2I5M13S', 163, 251, indels)

    def test_SMDMIM(self):
        indels = ([(6110, 1)], [(6120, 1)])
        self._test_read(6109, '66S2M1D9M1I173M', 83, 251, indels)

    def test_SMIMDM(self):
        indels = ([(7604, 1)], [(7605, 2)])
        self._test_read(7603, '12S2M1I1M2D235M', 83, 251, indels)

    def test_MDMIMDMS(self):
        indels = ([(202, 2)], [(199, 1), (210, 2)])
        self._test_read(1, '199M1D2M2I8M2D2M38S', 99, 251, indels)

    def test_SMDMDMDM(self):
        indels = ([], [(10530, 1), (10527, 1), (10522, 1)])
        self._test_read(10517, '32S6M1D4M1D2M1D207M', 147, 251, indels)

    def test_MIMIMDMS(self):
        indels = ([(13605, 2), (13614, 1)], [(13615, 2)])
        self._test_read(13388, '218M2I9M1I1M2D11M9S', 163, 251, indels)


def logging_setup():
    logging.basicConfig(stream=sys.stderr, level=logging.WARNING, format='%(message)s')
    for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
        level_name = logging.getLevelName(level)
        logging.addLevelName(level, level_name.capitalize())

if __name__ == '__main__':
  main(sys.argv)
