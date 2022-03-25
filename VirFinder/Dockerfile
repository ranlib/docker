FROM rocker/r-ubuntu:20.04
RUN apt-get update
RUN apt-get update --yes
RUN apt-get install --yes git
RUN git clone https://github.com/jessieren/VirFinder.git
RUN R -e "install.packages('glmnet', dependencies=TRUE)"
RUN R -e "install.packages('Rcpp', dependencies=TRUE)"
RUN R -e "install.packages('BiocManager')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(c('qvalue'))"
RUN R -e "install.packages('/VirFinder/linux/VirFinder_1.1.tar.gz', repos = NULL, type='source')"

