FROM ubuntu:24.04

RUN apt-get update && apt-get install -y \
    build-essential \
    python3 \
    python3-distutils-extra \
    python-is-python3 \
    git \
    zlib1g-dev

RUN git clone https://github.com/linsalrob/primer-trimming.git \
    && cd primer-trimming \
    && make all \
    && make install

#RUN git clone https://github.com/linsalrob/primer-trimming.git \
#    && cd primer-trimming \
#    && python3 setup.py install --user
