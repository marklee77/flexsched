FROM ubuntu:trusty
MAINTAINER Mark Stillwell <mark@stillwell.me>

ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && \
    apt-get -y install \
        build-essential \
        curl \
        cmake \
        git \
        libglpk-dev \
        pkg-config && \
    rm -rf /var/lib/apt/lists/* /var/cache/apt/*

COPY . /scratch/flexsched

WORKDIR /scratch
