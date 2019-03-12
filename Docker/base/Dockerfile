FROM jupyterhub/singleuser:0.9.4
MAINTAINER https://github.com/underworldcode/stripy
ENV IMAGENAME=underworldcode/py3-lavavu-notebook-bundle
ENV PROJ_NAME=py3-lavavu-notebook-bundle

## =============================================================
## base - image ... whatever functionality you want to provide !
## This is my unix / python stuff (but it doesn't have underworld)
##
## This dockerfile builds an image from this content, and serves the
## sample web pages and notebooks at port 8080
##
## docker run -p 8181:8080 --name="docker-web-notebooks-test" -t lmoresi/lmoresi/docker-web-notebooks-module
## and then browse the docker VM ip address on port 8181 (for example)
##
## OR just use kitematic and click on the preview image
##
## =============================================================

USER root

# The following installations are required for:
# various: build-essential
# stripy: gcc, gfortran
# lavavu: libgl*mesa*, libpng-dev, libtiff-dev, libx11-dev
# lavavu: freeglut*, xvfb


RUN apt-get update -y ;
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends --fix-missing \
    gcc \
    gfortran \
    libpng-dev \
    libtiff-dev \
    libx11-dev \
    build-essential \
    libgl1-mesa-dri \
    libgl1-mesa-glx \
    freeglut3 \
    freeglut3-dev \
    xvfb \
    && conda install -y numpy scipy cartopy gdal netcdf4 \
    && rm -rf /var/lib/apt/lists/*



RUN python3 -m pip install --no-cache-dir stripy lavavu

# RUN useradd --create-home --home-dir /home/jovyan --shell /bin/bash --user-group jovyan || true

# # Install Tini.. this is required because CMD (below) doesn't play nice with notebooks for some reason: https://github.com/ipython/ipython/issues/7062, https://github.com/jupyter/notebook/issues/334
# RUN curl -L https://github.com/krallin/tini/releases/download/v0.10.0/tini > tini && \
# echo "1361527f39190a7338a0b434bd8c88ff7233ce7b9a4876f3315c22fce7eca1b0 *tini" | sha256sum -c - && \
# mv tini /usr/local/bin/tini && \
# chmod +x /usr/local/bin/tini

# script for xvfb-run. all docker commands will effectively run under this via the entrypoint
RUN printf "#\041/bin/sh \n rm -f /tmp/.X99-lock && xvfb-run -s '-screen 0 1600x1200x16' \$@" >> /usr/local/bin/xvfbrun.sh && \
chmod +x /usr/local/bin/xvfbrun.sh


# Jovyan user / group
RUN addgroup jovyan || true
RUN usermod -a -G jovyan jovyan || true

WORKDIR /home/jovyan
USER jovyan


# expose notebook port and server port
EXPOSE 8888

WORKDIR /home/jovyan/


# note we use xvfb which to mimic the X display for lavavu
ENTRYPOINT ["/usr/local/bin/xvfbrun.sh"]