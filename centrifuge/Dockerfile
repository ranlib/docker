FROM ubuntu:22.04
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates wget make g++ \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN wget https://github.com/DaehwanKimLab/centrifuge/archive/refs/tags/v1.0.4.tar.gz && tar -xvf v1.0.4.tar.gz
WORKDIR /centrifuge-1.0.4
RUN make
RUN rm /v1.0.4.tar.gz
ENV PATH /centrifuge-1.0.4:$PATH
