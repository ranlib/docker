FROM ubuntu:latest
RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install -y wget
RUN wget --quiet http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz && tar -xvf hmmer-3.1b2-linux-intel-x86_64.tar.gz
ENV PATH ${PATH}:/hmmer-3.1b2-linux-intel-x86_64/binaries
RUN wget --quiet https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux -O /usr/local/bin/prodigal && chmod u+x /usr/local/bin/prodigal
RUN apt-get clean
RUN apt-get autoremove

