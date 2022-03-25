FROM ubuntu:latest
ENV DEBIAN_FRONTEND noninteractive
ENV LC_ALL en_US.UTF-8
ENV LC_CTYPE en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8
RUN apt-get update -y
RUN apt-get upgrade -y
RUN apt-get purge locales
RUN apt-get install -y locales
RUN locale-gen en_US.UTF-8
RUN dpkg-reconfigure locales
RUN apt-get install -y python3-pip python-is-python3 wkhtmltopdf
RUN pip3 install requests pandas openpyxl tabulate pdfkit jira
COPY plantqc.py /usr/local/bin
RUN apt-get clean
RUN apt-get autoremove

