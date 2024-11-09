FROM ubuntu:24.04

ARG DEBIAN_FRONTEND=noninteractive
ARG SAMTOOLS_VERSION=1.21
ARG BCFTOOLS_VERSION=1.21
ARG HTSLIB_VERSION=1.21

RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -y build-essential wget zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev

WORKDIR /usr/src

RUN wget --quiet https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	tar jxf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	rm samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
	cd samtools-${SAMTOOLS_VERSION} && \
	./configure --without-curses --prefix $(pwd) && \
	make
	
RUN wget --quiet https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    tar xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    rm -rf bcftools-${BCFTOOLS_VERSION}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VERSION} && \
    ./configure --without-curses --prefix $(pwd) && \
    make

RUN wget --quiet https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar xjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    rm -rf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure --without-curses --prefix $(pwd) && \
    make

RUN apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV PATH=${PATH}:/usr/src/samtools-${SAMTOOLS_VERSION}:/usr/src/bcftools-${BCFTOOLS_VERSION}:/usr/src/htslib-${HTSLIB_VERSION}
