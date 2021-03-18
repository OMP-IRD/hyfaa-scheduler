FROM continuumio/miniconda3 AS base

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

RUN ulimit -s unlimited

#install libraries
RUN apt-get -q update \
    && apt-get -q -y install build-essential gfortran cmake libnetcdf-dev libnetcdff-dev  \
    && apt-get -q clean \
    && rm -rf /var/lib/apt/lists/*

# installs common to builder & production
RUN conda install -y -c conda-forge python=3.7 numpy scipy netCDF4 pyyaml progress pandas pytest requests SALib ftputil

COPY install.sh install.sh
COPY README.md README.md
COPY src src
RUN ./install.sh /hyfaa
ENV PATH /hyfaa:/hyfaa/bin:${PATH}
ENV PYTHONPATH /hyfaa/lib/python3.7/site-packages:${PYTHONATH}

RUN mkdir /work
WORKDIR /work

CMD ["/bin/bash", "run.sh"]


