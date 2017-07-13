#!/usr/bin/env bash
dirname=$(dirname $0)
set -ue

echo -e "\tBAMRead.indel_at/cigar-varieties.bam:"
$dirname/unit-diff-tests.py BAMRead.indel_at $dirname/cigar-varieties.bam -i $dirname/cigar-varieties.indel_at.tsv | diff -s - $dirname/cigar-varieties.indel_at.out
