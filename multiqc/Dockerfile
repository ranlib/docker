FROM ubuntu:latest@sha256:626ffe58f6e7566e00254b638eb7e0f3b11d4da9675088f4781a50ae288f3322 
ENV DEBIAN_FRONTEND noninteractive 
RUN apt-get update
RUN apt-get install -y software-properties-common \
    wget \
    python3 \
    python3-dev \
    python3-pip \
    texlive-xetex \
    texlive-science \
    locales && \
    locale-gen en_US.UTF-8 && \
    apt-get clean && apt-get autoclean && rm -rf /var/lib/apt/lists/*
RUN apt-get install --reinstall texlive-xetex

RUN pip3 install "multiqc==1.11"
ENV VERSION=2.16.2
RUN wget -q https://github.com/jgm/pandoc/releases/download/$VERSION/pandoc-$VERSION-linux-amd64.tar.gz && tar -xvf pandoc-$VERSION-linux-amd64.tar.gz && rm pandoc-$VERSION-linux-amd64.tar.gz
RUN fmtutil-sys --all 
ENV LC_ALL='C.UTF-8' LANG='C.UTF-8'
ENV PATH="/pandoc-$VERSION/bin:${PATH}"
RUN mkdir /data
WORKDIR /data
