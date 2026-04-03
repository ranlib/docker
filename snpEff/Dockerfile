FROM ubuntu:25.04
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates wget unzip default-jre && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
#RUN wget -q https://snpeff.odsp.astrazeneca.com/versions/snpEff_latest_core.zip && unzip snpEff_latest_core.zip snpEff/snpEff.jar && rm snpEff_latest_core.zip
RUN wget -q https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip&& unzip snpEff_latest_core.zip snpEff/snpEff.jar && rm snpEff_latest_core.zip
RUN printf '#!/bin/bash \njava -jar /snpEff/snpEff.jar "$@"\n' > /usr/local/bin/snpEff && chmod ogu+x /usr/local/bin/snpEff

