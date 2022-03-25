FROM rocker/r-ubuntu:20.04
RUN apt-get update --yes
RUN apt-get install --yes build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
RUN apt-get install --yes r-bioc-biostrings r-cran-knitr
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'devtools::install_git("https://gitlab.com/sirarredondo/mlplasmids",repos="http://cran.us.r-project.org")'
ADD mlplasmids.R /usr/local/bin/mlplasmids.R
