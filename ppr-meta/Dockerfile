FROM ubuntu:bionic

#install dependencies for PPR-Meta-tool
#falls gpu vorhanden: pip install tensorflow-gpu==1.4.1  #GPU version
RUN apt-get update && apt install -y procps python2.7 python-pip git unzip wget libxt6 && \
    git clone https://github.com/zhenchengfang/PPR-Meta.git && \
    pip install -U numpy h5py tensorflow keras==2.0.8 && \
    apt-get clean &&  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

#installing MCR (MatLab Runtime Compiler) 
RUN mkdir /mcr-install && \
    mkdir /opt/mcr && \
    cd /mcr-install && \
    wget -q http://ssd.mathworks.com/supportfiles/downloads/R2018a/deployment_files/R2018a/installers/glnxa64/MCR_R2018a_glnxa64_installer.zip && \
    cd /mcr-install && \
    unzip -q MCR_R2018a_glnxa64_installer.zip && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf mcr-install

ENV LD_LIBRARY_PATH /opt/mcr/v94/runtime/glnxa64:/opt/mcr/v94/bin/glnxa64:/opt/mcr/v94/sys/os/glnxa64:/opt/mcr/v94/extern/bin/glnxa64

RUN chmod +x /PPR-Meta/PPR_Meta
ENV PATH /PPR-Meta:$PATH
ENV PYTHONPATH /PPR-Meta:$PYTHONPATH
ENV H5PY_DEFAULT_READONLY=1
