FROM ubuntu:22.04

LABEL author="daniel.carrera@bms.com, carlos.rios@bms.com"
LABEL description="Docker container containing Sparta tools from Lyubetskaya, Anna <Anna.Lyubetskaya@bms.com>"
LABEL version=1.1
LABEL repo=XXXX

ENV PATH /usr/bin:/sbin:/usr/sbin:/usr/local/bin:/bin:/usr/lib/x86_64-linux-gnu/pkgconfig
ENV PKG_CONFIG_PATH="${PKG_CONFIG_PATH}:/usr/lib/x86_64-linux-gnu/pkgconfig"
ENV JAVA_HOME /usr/local/jdk1.7.0_71
ENV R_VERSION=4.2.3
ENV PYTHON_VERSION=3.11.2
ARG DEBIAN_FRONTEND=noninteractive 

RUN apt upgrade -y
RUN apt-get update -y
# required to install R packages
RUN apt-get install -y gdebi-core wget curl libcurl4-openssl-dev libxml2-dev libfontconfig1-dev libssl-dev libharfbuzz-dev libfribidi-dev libudunits2-dev libjpeg-dev libtiff5-dev libpng-dev libfreetype6-dev libmagick++-dev ffmpeg libsm6 libxext6 
# required for spdep
RUN apt-get install -y libgdal-dev
# required for hdf5r
RUN apt-get install -y libhdf5-dev patch

# install R and required libraries
RUN apt install -y g++ g++-11 gfortran gfortran-11 libgfortran-11-dev libopenblas-dev libopenblas-pthread-dev libopenblas0 libopenblas0-pthread make unzip zip 
RUN curl -O https://cdn.rstudio.com/r/ubuntu-2204/pkgs/r-${R_VERSION}_1_amd64.deb
RUN gdebi -n ./r-${R_VERSION}_1_amd64.deb
ENV PATH="${PATH}:/opt/R/${R_VERSION}/bin"
RUN /opt/R/4.2.3/bin/R --version

# install R packages - per Anna, leave Seurat at latest
RUN Rscript -e "install.packages(c('Seurat', 'remotes', 'tidyverse','spatstat','reticulate','spdep','hdf5r','optparse','aws.s3', 'BiocManager'), repos=c('https://pm.rdcloud.bms.com/prod-cran/2023-10-03') )"
RUN Rscript -e "BiocManager::install(c('glmGamPoi'), dependencies=TRUE)"
RUN Rscript -e 'remotes::install_github("zijianni/SpotClean@8a531d3", repos = "https://pm.rdcloud.bms.com/prod-cran/2023-10-03")'

# install Python and packages
RUN curl -O https://cdn.rstudio.com/python/ubuntu-2204/pkgs/python-${PYTHON_VERSION}_1_amd64.deb
RUN gdebi -n ./python-${PYTHON_VERSION}_1_amd64.deb
ENV PATH="${PATH}:/opt/python/${PYTHON_VERSION}/bin"
RUN pip install --upgrade pip

COPY . /P02567_TBIO-3021_10X_pilot
WORKDIR /P02567_TBIO-3021_10X_pilot

RUN pip install --no-cache-dir -r requirements.txt