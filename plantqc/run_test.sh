#!/bin/bash
DOCKER="docker run --rm -v $PWD:/mnt -w /mnt dbest/plantqc:latest "

# pacbio
#plantqc.py -l HSZWT -g -e -o HSZWT.tsv -s pacbio -v 
#plantqc.py -l HSZWT -g -t -e -o HSZWT.tsv -s pacbio -v 
#diff HSZWT.tsv tests/plantqc.tsv
#${DOCKER}plantqc.py -l HSZWT -g -t -e -o HSZWT.tsv -s pacbio -v 

# illumina
#        plantqc.py -l HTNNP -g -t -e -o HTNNP.tsv -s illumina -v -w -m -p
${DOCKER}plantqc.py -l HTNNP -g -t -e -o HTNNP.tsv -s illumina -v -w -m -p


exit 0
