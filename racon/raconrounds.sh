#!/bin/bash
shopt -s expand_aliases

if [ $# -ne 5 ] ; then
    echo "Usage: $0 <pacbio ccs fastq.gz file> <reference genome> <output directory> <racon rounds> <number of threads>"
    exit 1
fi

FASTQ=$1
ASSEMBLY=$2
OUTPUT=${3:-out_racon}
ROUNDS=${4:-2}
THREADS=${5:-8}

if [[ $HOSTNAME == "cori"* ]] ; then
    DOCKER="time shifter -E -m none --image="
else
    DOCKER="time docker run --rm -v $PWD:/mnt -w /mnt "
fi

# test undocker
#time raconrounds.py -r $FASTQ -a $ASSEMBLY --output ${OUTPUT} --pacbio --rounds $ROUNDS  -t $THREADS --debug

# variant calling
alias samtools='docker run -it --rm -v `pwd`:/mnt -w /mnt dbest/samtools:1.13 samtools'

samtools view -@ 10 -m 10G -o ${OUTPUT}/round_1.bam ${OUTPUT}/round_1.sam
samtools sort -@ 10 -m 10G -o ${OUTPUT}/round_1.sorted.bam ${OUTPUT}/round_1.bam
samtools index  ${OUTPUT}/round_1.sorted.bam
time freebayes --min-coverage 10 -p 1 -b ${OUTPUT}/round_1.sorted.bam -f ${ASSEMBLY} -v ${OUTPUT}/round_1.vcf

samtools view -@ 10 -m 10G -o ${OUTPUT}/round_2.bam ${OUTPUT}/round_2.sam
samtools sort -@ 10 -m 10G -o ${OUTPUT}/round_2.sorted.bam ${OUTPUT}/round_2.bam
samtools index  ${OUTPUT}/round_2.sorted.bam
time freebayes --min-coverage 10 -p 1 -b ${OUTPUT}/round_2.sorted.bam -f ${OUTPUT}/round_1_polished.fasta -v ${OUTPUT}/round_2.vcf

count_snps.sh ${OUTPUT}/round_2.vcf

# test docker
#${DOCKER}dbest/raconrounds:latest raconrounds -t $THREADS -r $FASTQ -a $ASSEMBLY --pacbio --rounds $ROUNDS --output ${OUTPUT} 
#${DOCKER}dbest/raconrounds:latest raconrounds -t $THREADS -r $FASTQ -a $ASSEMBLY --pacbio --rounds $ROUNDS --output ${OUTPUT} --final
#${DOCKER}dbest/raconrounds:latest raconrounds -t $THREADS -r $FASTQ -a $ASSEMBLY --pacbio --rounds $ROUNDS --output ${OUTPUT} --debug

# turn on trimming
#${DOCKER}dbest/raconrounds:latest raconrounds --trimming -t $THREADS -r $FASTQ -a $ASSEMBLY --pacbio --rounds $ROUNDS --output ${OUTPUT} --debug

exit 0
