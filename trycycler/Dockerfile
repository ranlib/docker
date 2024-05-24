FROM condaforge/miniforge3
RUN mamba create -c bioconda -c conda-forge -n trycycler trycycler
ENV PATH=$PATH:/opt/conda/envs/trycycler/bin
RUN echo "source activate trycycler" > ~/.bashrc
