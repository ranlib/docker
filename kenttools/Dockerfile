FROM ubuntu:latest
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update
RUN apt-get install --yes rsync libkrb5-3 curl
WORKDIR /ucsc
RUN rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ./
ENV PATH $PATH:/ucsc
RUN apt-get clean && apt-get autoremove


