FROM continuumio/miniconda3:24.3.0-0
COPY ./environment.yaml /
RUN conda install -c conda-forge libgcc-ng
RUN conda env create -f /environment.yaml && conda clean -a
RUN conda install -c bioconda pypolca
