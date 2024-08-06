#!/bin/bash
#
# test fastp with some small public data
#
wget https://raw.githubusercontent.com/StaPH-B/docker-builds/master/tests/SARS-CoV-2/SRR13957123_1.fastq.gz
wget https://raw.githubusercontent.com/StaPH-B/docker-builds/master/tests/SARS-CoV-2/SRR13957123_2.fastq.gz

docker run --rm -u 1002:1112 -v $PWD:/mnt -w /mnt dbest/fastp:v0.23.4 fastp \
        -i SRR13957123_1.fastq.gz \
        -I SRR13957123_2.fastq.gz \
        -o SRR13957123_PE1.fastq.gz \
        -O SRR13957123_PE2.fastq.gz \
        -h SRR13957123_fastp.html \
        -j SRR13957123_fastp.json
