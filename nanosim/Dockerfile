FROM continuumio/miniconda3:24.3.0-0
RUN conda update conda \
    && conda install -c conda-forge mamba 
RUN conda init bash \
    && . ~/.bashrc \
    && conda create -n nanosim python=3.8 \
    && conda activate nanosim \
    && mamba install -c bioconda -c conda-forge scikit-learn=0.22.1 six samtools pysam pybedtools minimap2 joblib htseq genometools-genometools \
    && mamba install -c bioconda -c conda-forge nanosim=3.2.2 numpy=1.19.5
ENV PATH=/opt/conda/envs/nanosim/bin/:${PATH}
