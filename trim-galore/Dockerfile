FROM ubuntu:24.04
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update -y && apt-get install -y --no-install-recommends ca-certificates curl cutadapt fastqc 
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.10.tar.gz -o trim_galore.tar.gz && tar xvzf trim_galore.tar.gz && rm trim_galore.tar.gz
RUN apt-get autoclean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
ENV PATH /TrimGalore-0.6.10:$PATH
