FROM ubuntu:24.04
RUN apt-get update
RUN apt-get install -y –no-install-recommends build-essential wget libssl-dev zlib1g zlib1g-dev && rm -rf /var/lib/apt/lists/*
RUN wget --quiet https://github.com/Kitware/CMake/releases/download/v3.29.3/cmake-3.29.3.tar.gz && tar xf cmake-3.29.3.tar.gz && rm cmake-3.29.3.tar.gz
RUN cd /cmake-3.29.3 && ./configure && make && make install
RUN wget --quiet https://github.com/lbcb-sci/racon/archive/refs/tags/1.5.0.tar.gz && tar -xvf racon-v1.5.0.tar.gz && rm racon-v1.5.0.tar.gz
RUN cd /racon-v1.5.0 && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make
ENV PATH $PATH:/racon-v1.5.0/build/bin

