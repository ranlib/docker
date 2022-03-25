FROM ubuntu:latest
RUN apt-get update
RUN apt-get install -y git wget make g++ python3 python-is-python3
#RUN git clone https://github.com/infphilo/centrifuge
RUN wget https://github.com/DaehwanKimLab/centrifuge/archive/refs/tags/v1.0.4.tar.gz && tar -xvf v1.0.4.tar.gz
RUN cd /centrifuge-1.0.4 && make
ENV PATH /centrifuge-1.0.4:$PATH
