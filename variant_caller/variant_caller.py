#!/usr/bin/env python
"""
simple variant caller
Take most common base in pileup of reads at a certain position.
Compare to base in reference at that position.
"""
import os
import sys
import re
import argparse
import pysam
import collections
import csv

def bases(samfile, reference_index, contig_name, minbasequal=20, verbose=False):
    """ Return base calls for given contig """
    calls = []
    for pileupcolumn in samfile.pileup(reference=contig_name, quality_threshold=minbasequal):
        bases = ""
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                bases += pileupread.alignment.query_sequence[pileupread.query_position]

        frequencies = collections.Counter("".join(sorted(bases)))
        my_call = frequencies.most_common(1)
        if my_call:
            coverage = int(pileupcolumn.nsegments)
            position = int(pileupcolumn.reference_pos) # position in reference sequence, 0-based!
            ref_base = reference_index.fetch(contig_name, position, position+1)
            base_call = my_call[0][0]
            base_call_count = my_call[0][1]
            if base_call != ref_base:
                # make position 1-based!
                calls.append( ( contig_name, position+1, ref_base, base_call, base_call_count, coverage, bases )  )
                if verbose:
                    print(f"{contig_name}\t{position+1}\t{ref_base}\t{base_call}\t{base_call_count}\t{coverage}\t{bases}")
                    
    return calls

def get_chromosomes(bam, keep_contigs=None):
    """ Create one interval for each chromosome """
    if keep_contigs is None:
        keep_contigs = "."
    pattern = re.compile(keep_contigs)
    aln = pysam.AlignmentFile(bam, 'rb')
    idxstats = aln.get_index_statistics()
    keep_contigs = []
    for i in idxstats:
        if i.mapped > 0 and pattern.match(i.contig):
            keep_contigs.append(i.contig)
    conlen = {x: aln.get_reference_length(x) for x in keep_contigs}
    aln.close()
    return conlen

if __name__ == '__main__':
    " make variant calls "
    parser = argparse.ArgumentParser(description='snp_caller')
    parser.add_argument('-r', '--reference', dest='reference', type=str, help='reference sequence in fasta format', required=True)
    parser.add_argument('-a', '--alignment', dest='alignment', type=str, help='alignment file in bam format', required=True)
    parser.add_argument('-o', '--output', dest='output', type=str, help='output file of variant calls', required=True)
    parser.add_argument('-q', '--minbasequal', dest='minbasequal', type=int, help='minimum base quality (default=20)', default=20)
    parser.add_argument('-v', '--verbose', action='store_true', help='debugging output')
    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    # check input exists
    if not os.path.exists(args.alignment):
        sys.exit(f"<E> input alignment file does not exist.")
    if not os.path.exists(args.reference):
        sys.exit(f"<E> input reference file does not exist.")
    if not os.path.exists(args.reference + '.fai'):
        pysam.faidx(args.reference)
    # get fasta index, contigs, read bam file
    findex = pysam.FastaFile(args.reference)
    contigs = get_chromosomes(args.alignment)
    bamfile = pysam.AlignmentFile(args.alignment, 'rb')
    # loop over contigs and make calls
    calls = {}
    for contig in contigs:
        calls[contig] = bases(bamfile, findex, contig, args.minbasequal, args.verbose)
    bamfile.close()
    # write calls to csv file
    with open(args.output,'w') as out:
        csv_out = csv.writer(out, delimiter='\t')
        header = ('contig_name', 'position', 'ref_base', 'base_call', 'base_call_count', 'coverage', 'bases')
        csv_out.writerow(header)
        for contig in contigs:
            for call in calls[contig]:
                csv_out.writerow(call)

    
