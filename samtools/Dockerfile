FROM ubuntu:22.04

RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -y build-essential wget libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev
RUN apt-get clean
RUN apt-get purge
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /usr/src

RUN wget --quiet https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
	tar jxf samtools-1.18.tar.bz2 && \
	rm samtools-1.18.tar.bz2 && \
	cd samtools-1.18 && \
	./configure --prefix $(pwd) && \
	make

ENV PATH=${PATH}:/usr/src/samtools-1.18
