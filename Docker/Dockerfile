#################################################
#  Short docker file to distribute some notebooks
#################################################

ARG FROMIMG_ARG=brmather/py3-lavavu-notebook-bundle:0.9.2-alpine
FROM ${FROMIMG_ARG}

##################################################
# Non standard as the files come from the packages

USER root

### Stripy - Notebooks

ENV MODULE_DIR="stripy-src"
ADD --chown=jovyan:jovyan $MODULE_DIR $MODULE_DIR
RUN   cd $MODULE_DIR && python3 -m pip install --no-cache-dir --no-deps --upgrade .

### litho1pt0 - Notebooks

ENV MODULE_DIR="litho1pt0-src"
ADD --chown=jovyan:jovyan $MODULE_DIR $MODULE_DIR
RUN   cd $MODULE_DIR && python3 -m pip install --no-cache-dir --no-deps --upgrade .

# This file is too big to bundle with the pip installation
ADD  --chown=jovyan:jovyan $MODULE_DIR/litho1pt0/Notebooks/litho1pt0/Data  Notebooks/Data


RUN ipython3 -c 'import stripy; stripy.documentation.install_documentation(path="Notebooks")'
ADD  --chown=jovyan:jovyan Notebooks/0-StartHere.ipynb Notebooks/0-StartHere.ipynb


RUN ipython3 -c 'import litho1pt0; litho1pt0.documentation.install_documentation(path="Notebooks")'

# This file is too big to bundle with the pip installation
# ADD  --chown=jovyan:jovyan $MODULE_DIR/miller_alaskamoho_srl2018/Notebooks/ShadedRelief/GRAY_HR_SR_OB.tif Notebooks/ShadedRelief/GRAY_HR_SR_OB.tif


# change ownership of everything
ENV NB_USER jovyan
RUN chown -R jovyan:jovyan /home/jovyan
USER jovyan

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
## ADD --chown=jovyan:jovyan  $NB_DIR /home/jovyan/Notebooks

# Trust all notebooks
RUN find -name \*.ipynb  -print0 | xargs -0 jupyter trust

# expose notebook port server port
EXPOSE $NB_PORT

VOLUME /home/jovyan/$NB_DIR/user_data

# note we use xvfb which to mimic the X display for lavavu
ENTRYPOINT ["/sbin/tini", "--", "/usr/local/bin/xvfbrun.sh"]

# launch notebook
ADD --chown=jovyan:jovyan Docker/scripts/run-jupyter.sh scripts/run-jupyter.sh
CMD scripts/run-jupyter.sh
