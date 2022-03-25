FROM ubuntu:latest
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update
RUN apt-get upgrade --yes
RUN apt-get install --yes git python3 python3-pip cython zlib1g-dev
RUN pip3 install setuptools
RUN pip3 install -U cython
#RUN pip3 install -U scikit-learn
RUN pip3 install scikit-learn==0.19.1
RUN git clone https://github.com/Shamir-Lab/PlasClass.git && cd PlasClass && python3 setup.py install && chmod u+x /PlasClass/classify_fasta.py
ENV PATH $PATH:/PlasClass
RUN apt-get clean
RUN apt-get autoremove

