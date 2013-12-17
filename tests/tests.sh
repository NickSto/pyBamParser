#!/usr/bin/env bash
dirname=$(dirname $0)
set -ue
echo -e "\tBAMRead.get_indels/cigar-varieties.bam:"
$dirname/../lib/unit.test.py BAMRead.get_indels $dirname/cigar-varieties.bam | diff -s - $dirname/cigar-varieties.out
#../lib/runner.test.py cigar-varieties.bam | diff -s - cigar-varieties.out
