FROM ubuntu:24.04
RUN apt-get update
RUN apt-get install -y --no-install-recommends build-essential wget cmake git ca-certificates zlib1g zlib1g-dev
RUN rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN wget --quiet --no-check-certificate https://github.com/lbcb-sci/raven/archive/refs/tags/1.8.3.tar.gz && tar xvf 1.8.3.tar.gz && rm 1.8.3.tar.gz
WORKDIR /raven-1.8.3
RUN mkdir build && cmake -S ./ -B ./build -DRAVEN_BUILD_EXE=1 -DCMAKE_BUILD_TYPE=Release
RUN cmake --build build -j 4
ENV PATH=${PATH}:/raven-1.8.3/build/bin
