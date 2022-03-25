# Building from the latest version of perl Docker image.
# https://hub.docker.com/_/perl/
FROM perl:latest

LABEL maintainer "Davis Templeton" <davistempleton3@gmail.com>

# Defining the version of Circos to build with.
ENV CIRCOS_VERSION 0.69-6

# Updating and installing some packages, creating a user, and installing Perl Modules.
# Then downloading Circos.
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y --no-install-recommends --no-install-suggests \
        libgd2-xpm-dev \
        vim \
    && rm -rf /var/lib/apt/lists/*; \
    useradd \
      --home-dir /home/circos \
      --shell /bin/bash \
      --create-home \
      circos; \
    cpanm --no-wget --notest \
        Clone \
        Config::General \
        Font::TTF::Font \
        GD \
        GD::Polyline \
        List::MoreUtils \
        Math::Bezier \
        Math::Round \
        Math::VecStat \
        Params::Validate \
        Readonly \
        Regexp::Common \
        SVG \
        Set::IntSpan \
        Statistics::Basic \
        Text::Format; \
    wget \
        --quiet \
        --directory-prefix=/opt \
        http://circos.ca/distribution/circos-${CIRCOS_VERSION}.tgz \
    && tar -xzf \
        /opt/circos-${CIRCOS_VERSION}.tgz \
        -C /opt \
    && rm /opt/circos-${CIRCOS_VERSION}.tgz;

# Set the user the Docker image will use from here on out.
USER circos

# Set the default working directory.
WORKDIR /home/circos

# Add Circos bin directory to user's PATH.
ENV PATH /opt/circos-${CIRCOS_VERSION}/bin:$PATH

# Set a default command to use when starting a container.
CMD [ "/bin/bash" ]