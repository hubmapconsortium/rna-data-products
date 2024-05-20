FROM --platform=linux/amd64 rocker/tidyverse

ENV DEBIAN_FRONTEND=non-interactive

RUN apt-get update \
 && apt-get -y install \
    build-essential \
    libcairo2-dev \
    libcurl4-gnutls-dev \
    libncurses-dev \
    libgeos-dev \
    libglpk-dev \
    libgsl-dev \
    libmagick++-dev \
    libxml2-dev \
    libssl-dev \
    python3-dev \
    python3-pip \
 && rm -rf /var/lib/apt/lists/*

# Make sure a 'python' command is available so anndata will install
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10

WORKDIR /opt

COPY install_R_packages.R /opt

RUN Rscript /opt/install_R_packages.R \
 && rm /opt/install_R_packages.R

COPY bin/azimuth.R /opt

CMD ["/bin/bash"]
