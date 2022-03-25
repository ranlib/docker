FROM ubuntu:latest
RUN apt-get update
RUN apt-get upgrade
RUN apt-get install --yes git python3
RUN git clone https://github.com/rpetit3/assembly-scan.git
ENV PATH $PATH:/assembly-scan

