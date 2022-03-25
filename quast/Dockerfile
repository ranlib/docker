FROM ubuntu:latest
ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y build-essential wget zlib1g-dev pkg-config libfreetype6-dev libpng-dev perl openjdk-8-jdk python pip python-setuptools python-dev locales libgd-dev
RUN pip install matplotlib
RUN locale-gen en_US.UTF-8
RUN dpkg-reconfigure locales
RUN cpan -i App::cpanminus
RUN cpanm --no-wget --notest Clone Config::General Font::TTF::Font GD GD::Polyline List::MoreUtils Math::Bezier Math::Round Math::VecStat Params::Validate Readonly Regexp::Common SVG Set::IntSpan Statistics::Basic Text::Format
RUN wget --quiet https://github.com/ablab/quast/releases/download/quast_5.1.0rc1/quast-5.1.0rc1.tar.gz
RUN tar -xvf quast-5.1.0rc1.tar.gz && cd quast-5.1.0rc1 && ./install.sh
RUN cd /quast-5.1.0rc1/quast_libs/bwa && make
RUN cd /quast-5.1.0rc1/quast_libs/bedtools && make
ENV PATH $PATH:/quast-5.1.0rc1
WORKDIR /
RUN wget --quiet http://circos.ca/distribution/circos-0.69-9.tgz
RUN tar -xvf circos-0.69-9.tgz
ENV PATH $PATH:/circos-0.69-9/bin
RUN apt-get clean
RUN apt-get autoremove

