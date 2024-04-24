FROM ubuntu:22.04
ENV PILON_VERSION=1.24
RUN apt-get update
RUN apt-get install -y default-jre wget
RUN wget https://github.com/broadinstitute/pilon/releases/download/v${PILON_VERSION}/pilon-${PILON_VERSION}.jar
RUN echo "#!/bin/bash" > pilon && echo "exec java -jar /pilon-${PILON_VERSION}.jar """"$""@"""" " >> pilon && chmod +x pilon
ENV PATH ${PATH}:/
