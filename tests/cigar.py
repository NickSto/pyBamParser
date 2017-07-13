#!/usr/bin/env python
from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals
import re
import os
import sys
import errno
import logging
import argparse
# Path hack to make sure the local pyBamParser loads before the installed version.
script_dir = os.path.dirname(os.path.realpath(__file__))  # .
lib_path = os.path.join(os.path.dirname(script_dir), 'lib')  # ../lib
sys.path.insert(0, lib_path)
import pyBamParser.read

OP_CHARS_TO_INTS = {'M':0, 'I':1, 'D':2, 'N':3, 'S':4, 'H':5, 'P':6, '=':7, 'X':8}
OP_INTS_TO_CHARS = {0:'M', 1:'I', 2:'D', 3:'N', 4:'S', 5:'H', 6:'P', 7:'=', 8:'X'}

ARG_DEFAULTS = {'log':sys.stderr, 'volume':logging.ERROR}
DESCRIPTION = """"""


def make_argparser():

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.set_defaults(**ARG_DEFAULTS)

    parser.add_argument('coord', type=int)
    parser.add_argument('-c', '--cigar',
        help='The SAM CIGAR field for the alignment.')
    parser.add_argument('-p', '--pos', type=int,
        help='The SAM POS field for the alignment.')
    parser.add_argument('-l', '--length', type=int,
        help='The read\'s length.')
    parser.add_argument('-r', '--reverse', action='store_true',
        help='The read is on the reverse strand (the POS is the 5\' end of the read).')
    parser.add_argument('-R', '--read', action='store_true',
        help='Convert coord from ref to read coordinates.')
    parser.add_argument('-s', '--seq')
    parser.add_argument('-L', '--log', type=argparse.FileType('w'),
        help='Print log messages to this file instead of to stderr. Warning: Will overwrite the file.')
    parser.add_argument('-q', '--quiet', dest='volume', action='store_const', const=logging.CRITICAL)
    parser.add_argument('-v', '--verbose', dest='volume', action='store_const', const=logging.INFO)
    parser.add_argument('-D', '--debug', dest='volume', action='store_const', const=logging.DEBUG)

    return parser


def main(argv):

    parser = make_argparser()
    args = parser.parse_args(argv[1:])

    logging.basicConfig(stream=args.log, level=args.volume, format='%(message)s')
    tone_down_logger()

    if args.cigar and args.pos and args.length:
        cigar = args.cigar
        pos = args.pos
        seq = args.seq
        reverse = args.reverse
        length = args.length
    elif args.cigar or args.pos or args.length:
        fail('--cigar, --pos, and --length are all required, unless reading from stdin.')
    else:
        cigar, pos, seq, reverse, length = parse_sam_line(sys.stdin.readline())

    read = pyBamParser.read.BAMRead('12345678901234567890123456789012', None)

    cigar_list = split_cigar(cigar)
    blocks = read._get_contiguous_blocks(pos, cigar_list, reverse, length)

    logging.debug(blocks)

    if args.read:
        converted_coord = read._to_read_coord(blocks, args.coord)
        read_coord = converted_coord
    else:
        converted_coord = read._to_ref_coord(blocks, args.coord)
        read_coord = args.coord

    if seq:
        print(format_seq(seq, read_coord))
    print(converted_coord)


def parse_sam_line(line):
    if line == '':
        fail('Error: No input')
    fields = line.rstrip('\r\n').split('\t')
    flags = int(fields[1])
    length = len(fields[9])
    return fields[5], int(fields[3]), fields[9], flags & 16, length


def split_cigar(cigar):
    cigar_list = []
    bits = re.findall(r'\d+[A-Z=]', cigar)
    for bit in bits:
        op_char = bit[-1:]
        length = int(bit[:-1])
        op_int = OP_CHARS_TO_INTS[op_char]
        cigar_list.append((length, op_int))
    return cigar_list


def format_seq(seq, coord):
    return (seq[max(0,coord-10):coord-1]+'|'+
            seq[coord-1]+'|'+
            seq[coord:coord+10])


def tone_down_logger():
    """Change the logging level names from all-caps to capitalized lowercase.
    E.g. "WARNING" -> "Warning" (turn down the volume a bit in your log files)"""
    for level in (logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG):
        level_name = logging.getLevelName(level)
        logging.addLevelName(level, level_name.capitalize())


def fail(message):
    logging.critical(message)
    if __name__ == '__main__':
        sys.exit(1)
    else:
        raise Exception('Unrecoverable error')


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv))
    except IOError as ioe:
        if ioe.errno != errno.EPIPE:
            raise
