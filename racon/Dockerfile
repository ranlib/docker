FROM ubuntu:24.04
RUN apt-get update
RUN apt-get install -y --no-install-recommends build-essential wget libssl-dev zlib1g zlib1g-dev cmake git ca-certificates && rm -rf /var/lib/apt/lists/*
RUN wget --no-check-certificate --quiet https://github.com/lbcb-sci/racon/archive/refs/tags/1.5.0.tar.gz && tar -xf 1.5.0.tar.gz && rm 1.5.0.tar.gz
WORKDIR /racon-1.5.0
RUN mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make
ENV PATH $PATH:/racon-1.5.0/build/bin
