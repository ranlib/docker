FROM ubuntu:latest
RUN apt-get update
RUN apt-get install --yes wget
RUN wget https://github.com/shenwei356/seqkit/releases/download/v0.12.1/seqkit_linux_amd64.tar.gz -O /tmp/seqkit.tar.gz && tar zxvf /tmp/seqkit.tar.gz -C /usr/bin/ && rm /tmp/seqkit.tar.gz
