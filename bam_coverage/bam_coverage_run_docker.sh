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

#
# micro test
#

# unsorted bam test
# rm -fr test_input
# mkdir -p test_input
# cp test_data/contig_7.bam test_input
# ${DOCKER}dbest/coverage:latest bam_coverage.py -b test_input/contig_7.bam -o test_output_micro/docker/unsorted_bam -u 10 -l 5 -v -t $CPU --width 1 --gwidth 100 -y
# head test_output_micro/docker/unsorted_bam/coverage.tsv

# # sam test
# rm -fr test_input
# mkdir -p test_input
# cp test_data/contig_7.sam test_input
# ${DOCKER}dbest/coverage:latest bam_coverage.py -s test_input/contig_7.sam -o test_output_micro/docker/sam -u 10 -l 5 -v -t $CPU --width 1 --gwidth 100 -y 
# head test_output_micro/docker/sam/coverage.tsv 


# #
# # test
# #

# # unsorted bam test
# mkdir -p test_input
# cp test_data/aln.bam test_input
# ${DOCKER}dbest/coverage:latest coverage -b test_input/aln.bam -o test_output/docker/unsorted_bam -u 10 -l 5 -v -t $CPU --gwidth 100 -y
# head test_output/docker/unsorted_bam/coverage.tsv

# # sam test
# mkdir -p test_input
# cp test_data/aln.sam test_input
# ${DOCKER}dbest/coverage:latest coverage -s test_input/aln.sam -o test_output/docker/sam -u 10 -l 5 -v -t $CPU --gwidth 100 -y
# head test_output/docker/sam/coverage.tsv 

# unsorted bam test multi processing
mkdir -p test_input
cp test_data/aln.bam test_input
${DOCKER}dbest/coverage:latest coverage -b test_input/aln.bam -o test_output_multi/docker/unsorted_bam -u 10 -l 5 -v -t $CPU --gwidth 100 -p 8 -y
head test_output_multi/docker/unsorted_bam/coverage.tsv

exit 0

