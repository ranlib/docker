FROM ubuntu:latest
RUN apt-get update
RUN apt-get install --yes git wget python3 python3-pip
RUN git clone https://github.com/nv-tlabs/meta-sim.git
RUN python3 -m pip install pur
RUN cd /meta-sim && pur -r requirements.txt && python3 -m pip install -r requirements.txt
RUN cd /meta-sim && scripts/data/download_assets.sh
RUN cd /meta-sim && python3 scripts/data/generate_dataset.py --config data/generator/config/mnist_val.json
RUN cd /meta-sim && python3 scripts/data/generate_dataset.py --config data/generator/config/bigmnist_val.json
RUN apt-get clean && apt-get autoremove
RUN export PYTHONPATH=$PWD:$PYTHONPATH



