FROM ubuntu:latest
ENV TZ=Europe/Kiev
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update
RUN apt-get install --yes wget git cmake g++ zlib1g zlib1g-dev libbz2-1.0 libbz2-dev
RUN git clone --recursive https://github.com/martin-steinegger/conterminator 
RUN mkdir /conterminator/build && cd /conterminator/build && cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. .. && make -j 4 && make install
ENV PATH /conterminator/build/bin/:$PATH
