FROM ubuntu:latest@sha256:626ffe58f6e7566e00254b638eb7e0f3b11d4da9675088f4781a50ae288f3322
RUN apt-get update --yes
RUN apt-get upgrade --yes
RUN apt-get install --yes python3 python3-pip python-is-python3 git ncbi-blast+
RUN pip3 install biopython requests numpy pandas joblib scikit-learn==0.22.2.post1
#RUN git clone --recursive https://github.com/leaemiliepradier/PlasForest
COPY PlasForest /PlasForest
#COPY plasmid_refseq.fasta* /PlasForest
RUN cd /PlasForest && tar -zxvf plasforest.sav.tar.gz
RUN chmod a+rwx /PlasForest && chmod a+rwx /PlasForest/*.py
RUN chmod a+x /PlasForest/*.sh
ENV PATH $PATH:/PlasForest
RUn mkdir /db
