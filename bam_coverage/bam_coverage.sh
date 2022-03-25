#!/bin/bash
SMRTLINK=bryce911/smrtlink:10.1.0.119588
BIN=/smrtlink/install/smrtlink-release_10.1.0.119588/smrtcmds/bin

if [ $# -ne 2 ] ; then
    echo "Usage: $0 <assembly fasta file>  <pacbio ccs fastq.gz file>"
    exit 1
fi

OUTPUT=coverage

ASS=$1
FAQ=$2

if [[ $HOSTNAME == "cori"* ]] ; then
    DOCKER="time shifter -E --image="
    CPU=32
else
    DOCKER="time docker run --rm -v $PWD:/mnt -w /mnt "
    CPU=8
fi

# align with minimap2
mkdir -p ${OUTPUT}/minimap2_output
SAM=${OUTPUT}/minimap2_output/aln.sam
${DOCKER}staphb/minimap2 minimap2 -a -x map-pb -o $SAM $ASS $FAQ

# coverage
mkdir -p ${OUTPUT}/coverage_output
${DOCKER}dbest/coverage:v1.0 bam_coverage.py -s $SAM -o $OUTPUT/coverage_output -v

exit 0
