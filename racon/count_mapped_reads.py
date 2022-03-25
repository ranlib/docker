#!/usr/bin/env python

import sys
import pysam

bam_file = sys.argv[1]
mapped = full_length_mapped = 0
for i in pysam.AlignmentFile(bam_file, "r"):
    if i.is_unmapped or i.is_supplementary or i.is_secondary:
        continue
    qseq = i.query_sequence.upper()
    rseq = i.get_reference_sequence().upper()
    mapped += 1
    if qseq == rseq:
        full_length_mapped += 1

print('mapped: %d full_length_mapped: %d' % (mapped, full_length_mapped))

