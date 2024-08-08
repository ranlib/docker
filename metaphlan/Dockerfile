FROM continuumio/miniconda3:24.5.0-0
RUN conda install --yes -c bioconda biopython
RUN conda install --yes -c bioconda -c conda-forge metaphlan=4.1.1
#RUN conda install --yes --solver=classic -c default -c conda-forge unzip
#RUN wget --quiet https://github.com/BenLangmead/bowtie2/releases/download/v2.5.4/bowtie2-2.5.4-linux-x86_64.zip \
#    && unzip bowtie2-2.5.4-linux-x86_64.zip
#RUN rm bowtie2-2.5.4-linux-x86_64.zip
#RUN apt-get remove unzip
#RUN metaphlan --install
