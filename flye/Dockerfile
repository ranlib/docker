FROM ubuntu:24.04
ARG FLYE_VER="2.9.5"

RUN apt-get update && apt-get install -y make \
 gcc \
 g++ \
 python3 \
 python3-pip \
 python-is-python3 \	
 zlib1g-dev \
 wget && \
 rm -rf /var/lib/apt/lists/* && apt-get clean && apt-get autoclean

# install flye and make /data
RUN wget --quiet https://github.com/fenderglass/Flye/archive/${FLYE_VER}.tar.gz && \
  tar -xvf ${FLYE_VER}.tar.gz && \
  rm ${FLYE_VER}.tar.gz && \
  cd Flye-${FLYE_VER} && \
  python setup.py build && \
  python setup.py install && \
  mkdir /data

WORKDIR /data

# just in case for singularity compatibility
ENV LC_ALL=C