FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
apt-get install --yes --no-install-recommends wget git python3 python3-pip python-is-python3 && \
rm -rf /var/lib/apt/lists/*
RUN pip install numpy pyfaidx tqdm
RUN git clone https://github.com/mkpython3/Mutation-Simulator.git /opt/Mutation-Simulator
ENV PYTHONPATH=$PYTHONPATH:/opt/Mutation-Simulator
RUN printf '#!/bin/bash \
\npython -m mutation_simulator "$@"\n' > /usr/local/bin/mutation_simulator && chmod ogu+x /usr/local/bin/mutation_simulator
