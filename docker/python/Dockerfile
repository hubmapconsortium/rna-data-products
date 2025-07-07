FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

## System packages
RUN apt-get update \
 && apt-get -y install \
    autoconf \
    automake \
    bison \
    flex \
    git \
    libfreetype6-dev \
    libxml2-dev \
    libtool \
    m4 \
    pkg-config \
    python3 \
    python3-dev \
    python3-pip \
    zlib1g-dev \
    curl \
 && rm -rf /var/cache/apt/*

## Python packages
COPY docker/python/requirements.txt /opt/requirements.txt
RUN pip3 install -r /opt/requirements.txt \
 && rm -rf /root/.cache/pip

RUN mkdir /root/.ssh

COPY bin /opt
COPY data /opt/data

RUN chmod -R +x /opt

WORKDIR /opt