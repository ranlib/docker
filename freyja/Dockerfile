FROM ubuntu:22.04
ARG FREYJA_VERSION=2.0.3
ARG SAMTOOLS_VERSION=1.23.1
ARG HTSLIB_VERSION=1.23.1
ARG IVAR_VERSION=1.4.4
ARG USHER_VER=0.6.3

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y \
        python3.11 \
        python3.10-venv \
        python3.11-distutils \
        python3-pip && \
    apt-get install -y --no-install-recommends wget make gcc g++ ca-certificates zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev bzip2 autoconf automake \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN update-alternatives --install /usr/bin/python python /usr/bin/python3.11 1
RUN python --version

WORKDIR /opt/src

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip3 install PyYAML setuptools click pandas requests cvxpy matplotlib tqdm seaborn plotly pyarrow sphinx

RUN wget -q https://github.com/andersen-lab/Freyja/archive/refs/tags/v${FREYJA_VERSION}.tar.gz && \
    tar -xvf v${FREYJA_VERSION}.tar.gz && \
    rm v${FREYJA_VERSION}.tar.gz && \
    cd Freyja-${FREYJA_VERSION} && \
    make install
    
RUN wget --quiet https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    tar jxf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    rm samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VERSION} && \
    ./configure --without-curses --prefix $(pwd) && \
    make && make install
	
RUN wget --quiet https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar xjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    rm -rf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure --without-curses --prefix $(pwd) && \
    make && make install

RUN wget --quiet https://github.com/andersen-lab/ivar/archive/refs/tags/v${IVAR_VERSION}.tar.gz && \
    tar xf v${IVAR_VERSION}.tar.gz && \
    rm -rf v${IVAR_VERSION}.tar.gz && \
    cd ivar-${IVAR_VERSION} && \
    sed -i 's|$with_hts/include|$with_hts|' configure.ac && \
    sed -i 's|$with_hts/lib|$with_hts|' configure.ac && \	
    ./autogen.sh && \
    ./configure --with-hts=/opt/src/htslib-${HTSLIB_VERSION} && \
    make && make install

ENV PATH ${PATH}:/opt/src/samtools-${SAMTOOLS_VERSION}/bin:/opt/src/ivar-${IVAR_VERSION}/src
ENV LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/opt/src/htslib-${HTSLIB_VERSION}/lib
