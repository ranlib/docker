FROM ubuntu:latest
ENV DEBIAN_FRONTEND noninteractive
ENV fastp_version 0.23.1
RUN apt-get update -y && apt-get install -y \
    build-essential \
    libnss-sss \
    wget \
    unzip \
    zlib1g-dev \
    libisal-dev \
    libdeflate-dev

WORKDIR /usr/local/bin/
RUN wget --quiet https://github.com/OpenGene/fastp/archive/v${fastp_version}.zip
RUN unzip v${fastp_version}.zip
WORKDIR /usr/local/bin/fastp-${fastp_version}
RUN make
RUN make install
WORKDIR /usr/local/bin
