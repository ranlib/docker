FROM ubuntu:latest
RUN apt-get update
RUN apt-get install --yes build-essential git wget zlib1g-dev
RUN git clone https://github.com/vpc-ccg/haslr.git && cd haslr && make

