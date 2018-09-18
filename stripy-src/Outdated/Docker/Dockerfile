FROM lmoresi/unimelb-debian-python:1.05

## =============================================================
##
## Docker file to test / run stripy (including litho1pt0 that depends
## on stripy and can be used as a testing platform)
##
## =============================================================


## Needed for lavavu

RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
        xorg-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff-dev \
        libxft-dev \
        xvfb \
        freeglut3 \
        freeglut3-dev \
        libgl1-mesa-dri \
        libgl1-mesa-glx \
        xauth \
        python-tk \
        swig && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir /python_modules
WORKDIR   /python_modules

ADD src src
ADD f77-src f77-src
ADD stripy stripy
ADD MANIFEST.in MANIFEST.in
ADD setup.py setup.py

RUN python setup.py build && python setup.py install

## litho1pt0 is built entirely on stripy and is a useful learning tool / testing tool
## but the notebooks are in the source code so should be installed from github as
## above for stripy.
## RUN pip install litho1pt0


# Lavavu visualisation is used in the examples

RUN pip install --upgrade pip
RUN pip install lavavu

# script for xvfb-run.  all docker commands will effectively run under this via the entrypoint
RUN printf "#\041/bin/sh \n rm -f /tmp/.X99-lock && xvfb-run -a -s '-screen 0 1600x1200x16' \$@" >> /usr/local/bin/xvfbrun.sh && \
    chmod +x /usr/local/bin/xvfbrun.sh


# Create a non-privileged user to run the notebooks and switch to this user for the server

RUN useradd --create-home --home-dir /home/demonstration --shell /bin/bash --user-group demon
RUN chown -R demon:demon /home/demonstration

WORKDIR /home/demonstration

ADD --chown=demon:demon Docker/scripts/run-jupyter.py run-jupyter.py
ADD --chown=demon:demon Notebooks/SphericalTriangulations/*.ipynb Notebooks/SphericalTriangulations/
ADD --chown=demon:demon Notebooks/Data Notebooks/Data


# Persistent / Shared space outside the container
VOLUME /home/demonstration/external

EXPOSE 8080


# Now switch to the default user for the container

USER demon
ENV HOME=/home/demonstration
ENV SHELL=/bin/bash
ENV USER=demon

ENTRYPOINT ["/usr/local/bin/tini", "--", "xvfbrun.sh"]

CMD ./run-jupyter.py
