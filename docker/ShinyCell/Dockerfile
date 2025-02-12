FROM library/r-base:4.4.0
ENV DEBIAN_FRONTEND=non-interactive

#RUN dpkg -r --force-all libglib2.0-0

RUN apt-get update \
 && apt-get -y --allow-downgrades install \
    build-essential \
    libcurl4-gnutls-dev \
    curl \
    libncurses-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libglib2.0-0=2.81.2-1 \
    libglib2.0-0t64 \
    libgirepository-2.0-0 \
    libglib2.0-dev-bin-linux \
    libglib2.0-bin \
    libglib2.0-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgeos-dev \
    libglpk-dev \
    libgsl-dev \
    libhdf5-dev \
    libxml2-dev \
    python3 \
    python3-setuptools \
    python3-pip \
    python3-anndata \
    python3-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

RUN Rscript -e "install.packages('devtools',repos = 'http://cran.us.r-project.org')"

COPY install_R_packages.R /opt

RUN Rscript /opt/install_R_packages.R \
 && rm /opt/install_R_packages.R

COPY bin/make_shinycell.R /opt
COPY data/ensembl_to_symbol.json /opt

RUN chmod -R +x /opt

CMD ["/bin/bash"]