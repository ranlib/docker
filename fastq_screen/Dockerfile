FROM ubuntu:22.04
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && \
    apt-get install -y --no-install-recommends wget bwa ca-certificates libfindbin-libs-perl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

WORKDIR /opt
RUN wget --quiet https://github.com/StevenWingett/FastQ-Screen/archive/refs/tags/v0.15.3.tar.gz && \
	tar -xzvf v0.15.3.tar.gz && \
	rm v0.15.3.tar.gz
ENV PATH=$PATH:/opt/FastQ-Screen-0.15.3
