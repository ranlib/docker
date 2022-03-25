FROM ubuntu:latest
RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get install -y python3 python3-pip python-is-python3
RUN pip3 install pyhmmer==0.4.5
COPY AntiFam.* /AntiFam
COPY hmmer.py /usr/local/bin
RUN apt-get clean
RUN apt-get autoremove
