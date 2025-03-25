FROM ubuntu:24.04
ENV DEBIAN_FRONTEND="noninteractive"
RUN apt-get update -y
RUN apt-get install -y build-essential zlib1g-dev wget
RUN wget https://github.com/chhylp123/hifiasm/archive/refs/tags/0.25.0.tar.gz && tar -xvf 0.25.0.tar.gz && cd /hifiasm-0.25.0 && make
ENV PATH $PATH:/hifiasm-0.25.0
