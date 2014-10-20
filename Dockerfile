# Dockerizing Charite compbio NGS pipeline.
#
# Based on Ubuntu 14.04.

FROM       ubuntu:14.04
MAINTAINER M. Holtgrewe <manuel.holtgrewe@charite.de>

# Make sure to place the following line in ~/.profile (or execute on command line).
# This makes Docker use the Charite DNS server instead of Google's (8.8.8.8).
# Yay! for the Charite network setup.
#
# export DOCKER_OPTS="--dns 141.42.3.33 --dns 141.42.1.11"

# ---------------------------------------------------------------------------
# Installation
# ---------------------------------------------------------------------------

# Charite proxy configuration, Yay!
ENV http_proxy  http://proxy.charite.de:8080
ENV https_proxy http://proxy.charite.de:8080
ENV HTTP_PROXY  http://proxy.charite.de:8080
ENV HTTPS_PROXY http://proxy.charite.de:8080

# Edit /etc/apt/sources.list to use Charite apt-get services.
RUN perl -p -i -e 's/archive.ubuntu.com/ubuntu.charite.de/g' \
        /etc/apt/sources.list

# Update from Ubuntu repositories.
RUN apt-get update && apt-get upgrade -y
# HTML to PDF conversion and PDF juggling dependencies
RUN apt-get install -y pdfjam wkhtmltopdf
# Python dependencies for quality control script
RUN apt-get install -y python-jinja2
# Install dependencies.
RUN apt-get install -y wget unzip python-setuptools

# install the cbpipeline Python stuff
RUN \
    mkdir -p /tmp/cbpipeline && \
    cd /tmp/cbpipeline && \
    wget https://github.com/holtgrewe/cbqc/archive/v0.1.0.zip && \
    unzip v0.1.0.zip && \
    cd cbqc-0.1.0 && \
    python setup.py install && \
    rm -rf /tmp/cbpipeline
