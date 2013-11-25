#!/usr/bin/env bash
set -ue
echo -e "\tcigar-varieties.bam:"
../lib/runner.test.py cigar-varieties.bam | diff -s - cigar-varieties.out
