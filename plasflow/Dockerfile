FROM ubuntu:latest

RUN apt-get update && \
	ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && \
	apt-get install -y build-essential checkinstall python3 && \
	apt-get install -y ca-certificates software-properties-common wget libreadline-gplv2-dev && \
	apt-get install -y libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libgdbm-dev libc6-dev libbz2-dev && \
        apt-get install -y apt-transport-https default.jdk

ENV TZ=Europe/Lisbon

RUN cd /opt

#Install Anaconda

RUN wget https://repo.continuum.io/archive/Anaconda3-4.2.0-Linux-x86_64.sh
RUN bash Anaconda3-4.2.0-Linux-x86_64.sh -b -p /opt/anaconda
RUN rm Anaconda3-4.2.0-Linux-x86_64.sh
RUN /opt/anaconda/bin/conda update conda

#Install R 3.5.1

#RUN apt-get purge r-base* r-recommended r-cran-* && apt autoremove && apt update
#RUN /opt/anaconda/bin/conda install -c conda-forge r=3.5.1

RUN cd /bin

#Install tensorflow

RUN /opt/anaconda/bin/conda install -c conda-forge tensorflow
RUN /opt/anaconda/bin/conda uninstall protobuf
#RUN apt-get install -y python-pip 
#RUN pip install -U tensorflow

RUN cd /opt

#Install plasflow

RUN /opt/anaconda/bin/conda install -c bioconda plasflow

#Install intel-openmp 

RUN /opt/anaconda/bin/conda install -c anaconda intel-openmp 

ENV PATH="/opt/anaconda/bin/:$PATH"

WORKDIR /data
