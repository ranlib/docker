FROM ubuntu:24.04
SHELL ["/bin/bash", "-c"]
RUN apt-get update && \
    apt-get install -y python3 python3-venv && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN python3 -m venv ~/sniffles2 \
    && source ~/sniffles2/bin/activate \
    && python3 -m pip install sniffles pysam
RUN echo "source /root/sniffles2/bin/activate" > ~/.bashrc
ENV PATH $PATH:/root/sniffles2/bin