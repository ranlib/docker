FROM ubuntu:22.04
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates build-essential zlib1g-dev wget && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN wget -q https://github.com/DMU-lilab/pTrimmer/archive/refs/tags/V1.3.4.tar.gz \
    && tar -xvf V1.3.4.tar.gz \
    && cd pTrimmer-1.3.4 \
    && make \
    && ln -s pTrimmer-1.3.4 ptrimmer \
    && cd .. \
    && rm V1.3.4.tar.gz
ENV PATH $PATH:/pTrimmer-1.3.4
