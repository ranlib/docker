FROM ubuntu:24.04
ARG mm2_version=2.30
ARG samtools_version=1.22.1

# install packages
RUN apt-get update --allow-releaseinfo-change
RUN apt-get install -y --no-install-recommends python3 wget bzip2 git ssh ca-certificates autoconf make cmake gcc g++ zlib1g-dev libcurl4-openssl-dev liblzma-dev libbz2-dev libdeflate-dev libssl-dev curl git-lfs \
    && rm -rf /var/lib/apt/lists/*

# copy other resources
COPY prepare_run.py /usr/local/bin/
COPY merge_read_group_tags.py /usr/local/bin/
COPY cat_as_fastq.py /usr/local/bin/

# install minimap2
RUN wget --quiet https://github.com/lh3/minimap2/releases/download/v${mm2_version}/minimap2-${mm2_version}_x64-linux.tar.bz2 \
        && tar xjf minimap2-${mm2_version}_x64-linux.tar.bz2 \
        && rm minimap2-${mm2_version}_x64-linux.tar.bz2 \
        && cd minimap2-${mm2_version}_x64-linux \
        && cp minimap2 /usr/local/bin/

# install samtools
RUN wget --quiet https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make \
    && make install

