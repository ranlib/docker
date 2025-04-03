FROM ubuntu:24.04
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y wget python3 python-is-python3 python3-pip 
RUN pip install --break-system-packages paralleltask
RUN wget --quiet https://github.com/Nextomics/NextDenovo/archive/refs/tags/2.5.2.tar.gz && tar -vxzf 2.5.2.tar.gz && rm 2.5.2.tar.gz
ENV PATH $PATH:/NextDenovo-2.5.2
RUN apt-get clean
RUN apt-get autoremove

