FROM continuumio/miniconda
RUN pip install numpy networkx nose pysam
RUN git clone https://github.com/Shamir-Lab/Recycler.git
RUN cd Recycler && git checkout tags/v0.62 && python setup.py install --user
RUN apt-get autoremove
RUN apt-get autoclean
ENV PATH $PATH:/Recycler/bin


