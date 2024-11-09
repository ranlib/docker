FROM ubuntu:24.04
RUN apt-get update
RUN apt-get install --yes --no-install-recommends ca-certificates wget && rm -rf /var/lib/apt/lists/*
RUN wget --quiet https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz -O /tmp/seqkit.tar.gz && tar -xvf /tmp/seqkit.tar.gz -C /usr/bin && rm /tmp/seqkit.tar.gz
