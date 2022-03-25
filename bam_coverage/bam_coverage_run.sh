#!/bin/bash

if [[ $HOSTNAME == "cori"* ]] ; then
    DOCKER="time shifter -E -m none --image="
    CPU=32
else
    DOCKER="time docker run --rm -v $PWD:/mnt -w /mnt "
    CPU=8
fi

#
# cleanup
#
rm -fr test_input
rm -fr test_output

# print help
bam_coverage.py --help

# reference output
head -20 testoutput/coverage.yaml

#
# micro test
#

# unsorted bam test
rm -fr test_input
mkdir -p test_input
cp test_data/contig_7.bam test_input
time bam_coverage.py --bam test_input/contig_7.bam -o test_output_micro/plain/unsorted_bam --upper 10 --lower 5 -v --threads 8 -y
head test_output_micro/plain/unsorted_bam/coverage.tsv

# sam test
rm -fr test_input
mkdir -p test_input
cp test_data/contig_7.sam test_input
time bam_coverage.py --sam test_input/contig_7.sam -o test_output_micro/plain/sam --upper 10 --lower 5 -v --threads 8 -y
head test_output_micro/plain/sam/coverage.tsv

exit

# build docker image latest
make -f bam_coverage.makefile build

# unsorted bam test
rm -fr test_input
mkdir -p test_input
cp test_data/contig_7.bam test_input
${DOCKER}dbest/coverage:latest bam_coverage.py -b test_input/contig_7.bam -o test_output_micro/docker/unsorted_bam -u 10 -l 5 -v -t $CPU --width 1 --gwidth 100 -y
head test_output_micro/docker/unsorted_bam/coverage.tsv

# sam test
rm -fr test_input
mkdir -p test_input
cp test_data/contig_7.sam test_input
${DOCKER}dbest/coverage:latest bam_coverage.py -s test_input/contig_7.sam -o test_output_micro/docker/sam -u 10 -l 5 -v -t $CPU --width 1 --gwidth 100 -y 
head test_output_micro/docker/sam/coverage.tsv 


#
# test 
#

# unsorted bam test
mkdir -p test_input
cp test_data/aln.bam test_input
time bam_coverage.py -b test_input/aln.bam -o test_output/plain/unsorted_bam --upper 10 --lower 5 -v --threads $CPU --gwidth 100 -y
head test_output/plain/unsorted_bam/coverage.tsv

# multi processing
mkdir -p test_input
cp test_data/aln.bam test_input
time bam_coverage.py -b test_input/aln.bam -o test_output_multi/plain/unsorted_bam --upper 10 --lower 5 -v --threads 8 --gwidth 100 -p 8 -y
head test_output_multi/plain/unsorted_bam/coverage.tsv

exit

# sam test
mkdir -p test_input
cp test_data/aln.sam test_input
time bam_coverage.py -s test_input/aln.sam -o test_output/plain/sam --upper 10 --lower 5 -v --threads $CPU --gwidth 100 -y
head test_output/plain/sam/coverage.tsv

#
# docker test
#

# build docker image
make -f bam_coverage.makefile build

# unsorted bam test
mkdir -p test_input
cp test_data/aln.bam test_input
${DOCKER}dbest/coverage:latest coverage -b test_input/aln.bam -o test_output/docker/unsorted_bam -u 10 -l 5 -v -t $CPU --gwidth 100 -y
head test_output/docker/unsorted_bam/coverage.tsv

# multi processing
mkdir -p test_input
cp test_data/aln.bam test_input
${DOCKER}dbest/coverage:latest coverage -b test_input/aln.bam -o test_output_multi/docker/unsorted_bam -u 10 -l 5 -v -t $CPU --gwidth 100 -p 8 -y
head test_output_multi/docker/unsorted_bam/coverage.tsv

# sam test
mkdir -p test_input
cp test_data/aln.sam test_input
${DOCKER}dbest/coverage:latest coverage -s test_input/aln.sam -o test_output/docker/sam -u 10 -l 5 -v -t $CPU --gwidth 100 -y
head test_output/docker/sam/coverage.tsv 

#
# compare output
#
echo "==> start compare output with diff"
diff test_output/docker/unsorted_bam/coverage.tsv test_output/docker/sam/coverage.tsv
diff test_output/plain/unsorted_bam/coverage.tsv  test_output/plain/sam/coverage.tsv
echo "==> end compare output with diff"

#
# cleanup
#
#rm -r test_input
#rm -r test_output

#
# check code
#
#pylint bam_coverage.py
#python -m py_compile bam_coverage.py

#
# format python code
#
#black -v --diff --line-length 1000 bam_coverage.py
#black -v --line-length 1000 bam_coverage.py

exit 0

