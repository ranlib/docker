FROM replikation/debian_basic
MAINTAINER Christian Brandt <christian.jena@gmail.com>
# install conda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN bash Miniconda3-latest-Linux-x86_64.sh -b -p /conda
RUN rm -f Miniconda3-latest-Linux-x86_64.sh
ENV PATH /conda/bin:$PATH

RUN conda install python=3.6 numpy theano keras
RUN git clone https://github.com/jessieren/DeepVirFinder
RUN chmod +x /DeepVirFinder/dvf.py
ENV PATH /DeepVirFinder:$PATH
ENTRYPOINT ["dvf.py"]
CMD ["-h"]
