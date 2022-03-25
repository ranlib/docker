FROM ubuntu:latest@sha256:adf73ca014822ad8237623d388cedf4d5346aa72c270c5acc01431cc93e18e2d
RUN apt-get update
RUN apt-get install --yes python3-pip python-is-python3
RUN pip3 install matplotlib seaborn pysam pyyaml
COPY ./bam_coverage.py /usr/local/bin
RUN ln -s /usr/local/bin/bam_coverage.py /usr/local/bin/coverage
RUN apt-get clean && apt-get autoremove

