FROM ubuntu:24.04

RUN apt-get update && apt-get install -y \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    wget

WORKDIR /opt/pear
RUN wget --quiet https://github.com/ressy/pear/archive/refs/tags/v0.9.11.tar.gz \
    && tar -xvf v0.9.11.tar.gz \
    && rm v0.9.11.tar.gz \	
    && cd pear-0.9.11 \
    && ./configure \
    && make \
    && make install

ENV PATH=/usr/local/bin:${PATH}

