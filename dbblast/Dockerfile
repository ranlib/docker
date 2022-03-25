FROM ubuntu:latest
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update
RUN apt-get upgrade --yes
RUN apt-get install --yes --no-install-recommends python3 python3-biopython python-is-python3 ncbi-blast+
COPY ./blast.py /usr/local/bin/dbblast
RUN apt-get clean
RUN apt-get autoremove
