FROM ubuntu:24.04
RUN apt-get update -y
RUN apt-get install -y wget unzip
RUN wget --quiet https://github.com/wdecoster/chopper/releases/download/v0.8.0/chopper-linux.zip \
    && unzip chopper-linux.zip \
    && chmod u+x chopper \
    && mv chopper /usr/local/bin/
