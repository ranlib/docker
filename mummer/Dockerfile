FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get -y upgrade && \
        apt-get install -y build-essential git wget python pip\
	autoconf automake gcc libtool \
	yaggo fig2dev gnuplot xfig

WORKDIR /usr/src

RUN git clone --recursive https://github.com/jgurtowski/ectools.git
RUN wget --quiet https://github.com/mummer4/mummer/archive/v4.0.0rc1.tar.gz && \
	tar xvfz v4.0.0rc1.tar.gz && \
	cd mummer-4.0.0rc1/ && \
	autoreconf -fi && \
	./configure --prefix=/usr/bin && \
	make && make install && \
	cd .. && rm v4.0.0rc1.tar.gz

ENV PATH="/usr/src/mummer-4.0.0rc1:${PATH}"

RUN chmod a+x /usr/src/ectools/*[a-z].py && \
	ln -s /usr/src/ectools/*[a-z].py /usr/bin/

VOLUME /home/data
WORKDIR /home/data

RUN pip install mummer-idotplot

RUN apt-get clean
RUN apt-get autoremove
