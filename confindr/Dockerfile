FROM ubuntu:latest
RUN apt-get update
RUN apt-get install -y wget
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-py37_4.10.3-Linux-x86_64.sh
RUN bash Miniconda3-py37_4.10.3-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-py37_4.10.3-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda update -y conda
RUN conda install -c bioconda confindr
#RUN conda install -c bioconda confindr
#RUN apt-get autoremove
#RUN apt-get autoclean


