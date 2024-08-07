FROM ubuntu:24.04
RUN apt-get update -y
RUN apt-get install -y wget python3 openjdk-8-jre-headless && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get autoclean

RUN wget --quiet https://sourceforge.net/projects/bbmap/files/BBMap_39.08.tar.gz && \
    tar -xzf BBMap_39.08.tar.gz && \
    rm BBMap_39.08.tar.gz

ENV PATH="${PATH}:/bbmap" LC_ALL=C
    
RUN wget --quiet https://github.com/ncbi/sra-human-scrubber/releases/download/2.2.1/sra-human-scrubber-2.2.1.tar.gz \
    && tar -xvf sra-human-scrubber-2.2.1.tar.gz \
    && rm sra-human-scrubber-2.2.1.tar.gz
    
ENV PATH="/usr/bin:/sra-human-scrubber-2.2.1/bin:/sra-human-scrubber-2.2.1/scripts:${PATH}"

WORKDIR /sra-human-scrubber-2.2.1

RUN ./init_db.sh
RUN ldconfig
RUN scrub.sh -p 6 -t
