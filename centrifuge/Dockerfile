FROM ubuntu:22.04
ARG VERSION=1.0.4.1
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates wget make g++ \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN wget --quiet https://github.com/DaehwanKimLab/centrifuge/archive/refs/tags/v${VERSION}.tar.gz && tar -xvf v${VERSION}.tar.gz
WORKDIR /centrifuge-${VERSION}
RUN make
RUN rm /v${VERSION}.tar.gz
ENV PATH /centrifuge-${VERSION}:$PATH
