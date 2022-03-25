FROM ubuntu:latest@sha256:9d6a8699fb5c9c39cf08a0871bd6219f0400981c570894cd8cbea30d3424a31f
RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install -y build-essential wget libssl-dev zlib1g zlib1g-dev python-is-python3 python3-pip
RUN pip3 install pysam

ARG CMAKE_VER="3.21.2"
RUN wget --quiet http://www.cmake.org/files/v3.21/cmake-${CMAKE_VER}-linux-x86_64.tar.gz && tar -xvf cmake-${CMAKE_VER}-linux-x86_64.tar.gz
ENV PATH $PATH:/cmake-${CMAKE_VER}-linux-x86_64/bin

ARG RACON_VER="1.4.21"
RUN wget --quiet https://github.com/lbcb-sci/racon/releases/download/${RACON_VER}/racon-v${RACON_VER}.tar.gz && tar -xvf racon-v${RACON_VER}.tar.gz 
RUN cd /racon-v${RACON_VER} && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make
ENV PATH $PATH:/racon-v${RACON_VER}/build/bin

ARG MINIMAP2_VER="2.22"
RUN apt-get install -y python bzip2 
RUN wget --quiet https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VER}/minimap2-${MINIMAP2_VER}_x64-linux.tar.bz2 && tar -jxvf minimap2-${MINIMAP2_VER}_x64-linux.tar.bz2
ENV PATH "${PATH}:/minimap2-${MINIMAP2_VER}_x64-linux"

COPY ./raconrounds.py /usr/local/bin
RUN ln -s /usr/local/bin/raconrounds.py /usr/local/bin/raconrounds

RUN apt-get autoclean && rm -rf /var/lib/apt/lists/*



