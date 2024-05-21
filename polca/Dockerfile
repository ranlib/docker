FROM ubuntu:24.04
RUN apt-get update --yes
RUN apt-get install --yes --no-install-recommends build-essential wget python3 pipx samtools bwa
RUN wget --quiet https://github.com/freebayes/freebayes/archive/refs/tags/v1.3.6.tar.gz && tar -xvf v1.3.6.tar.gz && rm v1.3.6.tar.gz
#WORKDIR /freebayes-1.3.6
#RUN make
RUN wget --quiet https://github.com/gbouras13/pypolca/archive/refs/tags/v0.3.1.tar.gz && tar -xvf v0.3.1.tar.gz && rm v0.3.1.tar.gz
WORKDIR /pypolca-0.3.1
RUN pipx install -e .
ENV PATH $PATH:/root/.local/bin
