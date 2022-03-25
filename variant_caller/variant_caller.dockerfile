FROM ubuntu:latest
RUN apt-get update
RUN apt-get install --yes python3-pip python-is-python3
RUN pip3 install matplotlib seaborn pysam
COPY ./variant_caller.py /usr/local/bin
RUN apt-get clean
RUN apt-get autoremove

