#################################################
#  Short docker file to distribute some notebooks
#################################################


# Note - we start from the standard
# stripy docker and re-build with the
# contents of the repository 

ARG FROMIMG_ARG=underworldcode/stripy:2.0
FROM ${FROMIMG_ARG}

################################################## 
# Non standard as the files come from the packages
##################################################

USER root

### Stripy - Notebooks

ENV MODULE_DIR="stripy-src"
ADD --chown=jovyan:jovyan . / stripy-src/
RUN   cd $MODULE_DIR && python3 -m pip install --no-cache-dir --no-deps --upgrade .


# This should overwrite any existing Notebooks directory

RUN ipython3 -c 'import stripy; stripy.documentation.install_documentation(path="Notebooks")'
RUN ipython3 -c 'import litho1pt0; litho1pt0.documentation.install_documentation(path="Notebooks")'

# This file is too big to bundle with the pip installation

RUN mkdir -p /usr/local/files && chown -R jovyan:jovyan /usr/local/files
ADD --chown=jovyan:jovyan Docker/scripts  /usr/local/files
ENV PATH=/usr/local/files:${PATH}


# change ownership of everything
ENV NB_USER jovyan
RUN chown -R jovyan:jovyan /home/jovyan
USER jovyan


# Run PyTest
############

RUN cd $MODULE_DIR/stripy/tests && pytest *.py

# Non standard as the files come from the packages
##################################################


## These are supplied by the build script
## build-dockerfile.sh

ARG IMAGENAME_ARG
ARG PROJ_NAME_ARG=stripy
ARG NB_PORT_ARG=8888
ARG NB_PASSWD_ARG=""
ARG NB_DIR_ARG="Notebooks"
ARG START_NB_ARG="0-StartHere.ipynb"

# The args need to go into the environment so they
# can be picked up by commands/templates (defined previously)
# when the container runs

ENV IMAGENAME=$IMAGENAME_ARG
ENV PROJ_NAME=$PROJ_NAME_ARG
ENV NB_PORT=$NB_PORT_ARG
ENV NB_PASSWD=$NB_PASSWD_ARG
ENV NB_DIR=$NB_DIR_ARG
ENV START_NB=$START_NB_ARG


## NOW INSTALL NOTEBOOKS

# (This is not standard - nothing to do here )

## The notebooks (and other files we are serving up)


# Trust all notebooks
RUN find -name \*.ipynb  -print0 | xargs -0 jupyter trust

# expose notebook port server port
EXPOSE $NB_PORT

VOLUME /home/jovyan/$NB_DIR/user_data

# note we use xvfb which to mimic the X display for lavavu
ENTRYPOINT ["/usr/local/bin/xvfbrun.sh"]

# launch notebook
ADD --chown=jovyan:jovyan Docker/scripts/run-jupyter.sh scripts/run-jupyter.sh
CMD scripts/run-jupyter.sh
