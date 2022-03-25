FROM ubuntu:20.04
RUN apt-get update
RUN apt-get install -y build-essential git wget libssl-dev zlib1g
RUN wget --quiet http://www.cmake.org/files/v3.20/cmake-3.21.2.tar.gz && tar xf cmake-3.21.2.tar.gz
RUN cd /cmake-3.21.2 && ./configure && make && make install
RUN apt-get install zlib1g-dev
RUN wget --quiet https://github.com/lbcb-sci/racon/releases/download/1.4.21/racon-v1.4.21.tar.gz && tar -xvf racon-v1.4.21.tar.gz 
RUN cd /racon-v1.4.21 && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make
ENV PATH $PATH:/racon-v1.4.21/build/bin

