FROM ubuntu:latest

RUN apt-get update && apt-get -y upgrade && apt-get install -y build-essential wget zlib1g-dev libbz2-dev liblzma-dev git && \
        apt-get clean && apt-get purge && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /usr/src

ARG CMAKE_VER="3.21.2"
RUN wget --quiet http://www.cmake.org/files/v3.21/cmake-${CMAKE_VER}-linux-x86_64.tar.gz && tar -xvf cmake-${CMAKE_VER}-linux-x86_64.tar.gz
ENV PATH $PATH:/usr/src/cmake-${CMAKE_VER}-linux-x86_64/bin

RUN wget --quiet https://github.com/vcflib/vcflib/releases/download/v1.0.2/vcflib-1.0.2-src.tar.gz && \
    tar -xvf vcflib-1.0.2-src.tar.gz && \
    cd vcflib &&  \
    mkdir -p build && \
    cd build && \
    cmake .. && \
    cmake --build . && \
    cmake --install .
    
RUN git clone --recursive https://github.com/ekg/freebayes.git && \
	cd freebayes && \
	git checkout tags/v9.9.13 && \
	make && \
	make install
