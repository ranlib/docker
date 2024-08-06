FROM ubuntu:24.04
ENV DEBIAN_FRONTEND noninteractive
ENV fastp_version 0.23.4
RUN apt-get update -y && apt-get install -y --no-install-recommends wget ca-certificates
RUN mkdir fastp && cd fastp && wget --quiet http://opengene.org/fastp/fastp.0.23.4 && mv fastp.0.23.4 fastp && chmod a+x ./fastp
RUN apt-get autoclean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
ENV PATH="/fastp:$PATH"