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
    libhdf5-dev \
    libmagick++-dev \
    libxml2-dev \
    libssl-dev \
    python3-dev \
    python3-pip \
 && rm -rf /var/lib/apt/lists/*

# Make sure a 'python' command is available so anndata will install
RUN update-alternatives --install /usr/bin/python python /usr/bin/python3 10

RUN pip install anndata==0.7.6

WORKDIR /opt

COPY install_R_packages.R /opt

RUN Rscript /opt/install_R_packages.R \
 && rm /opt/install_R_packages.R

COPY bin/make_shinycell.R /opt
RUN chmod -R +x /opt

CMD ["/bin/bash"]