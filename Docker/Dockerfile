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
WORKDIR /python_modules

ADD src src
ADD stripy stripy
ADD MANIFEST.in MANIFEST.in
ADD setup.py setup.py

RUN python setup.py build && python setup.py install

## This is built entirely on stripy and is a useful learning tool / testing tool
RUN pip install litho1pt0 
RUN pip install lavavu

## These are the build templates etc
# ADD Data/Reference Data/Reference
# ADD Data/Modules python_modules
# ADD mkdocs.yml mkdocs.yml
# ADD jupyter-server-theme jupyter-server-theme
# ADD scripts scripts
# ADD Documentation Documentation
# ADD Notebooks Notebooks
# ADD docs docs

# Stripy / Litho1pt0

# RUN cd python_modules && \
#        tar -xzf litho1pt0.tar.gz && \
#        tar -xzf stripy.tar.gz && \
#        cd ..

# RUN cd python_modules/litho1pt0 && \
#        python setup.py build && \
#        python setup.py install && \
#        cd ../..

# RUN cd python_modules/stripy && \
#       python setup.py build && \
#       python setup.py install && \
#       cd ../..

## Update / build the site

# Create a non-privileged user to run the notebooks and switch to this user for the server

RUN useradd --create-home --home-dir /home/demonstration --shell /bin/bash --user-group demon
RUN chown -R demon:demon /home/demonstration

WORKDIR /home/demonstration


USER demon
ENV HOME=/home/demonstration
ENV SHELL=/bin/bash
ENV USER=demon



# Launch the notebook server from the Notebook directory
# but perhaps there is something else that would do this.

# Persistent / Shared space outside the container

VOLUME /home/demonstration/external

ADD Docker/scripts/run-jupyter.py run-jupyter.py
ADD Notebooks Notebooks

EXPOSE 8080
ENTRYPOINT ["/usr/local/bin/tini", "--"]


CMD ./run-jupyter.py


# # expose notebook port
# EXPOSE 8888
# # expose glucifer port
# EXPOSE 9999

# # add default user jovyan and change permissions on NB_WORK
# ENV NB_USER jovyan
# RUN useradd -m -s /bin/bash -N jovyan

# # copy this file over so that no password is required
# COPY jupyter_notebook_config.json /home/$NB_USER/.jupyter/jupyter_notebook_config.json

# # update all permissions for user
# RUN chown -R $NB_USER:users $NB_WORK $UW2_DIR /home/$NB_USER

# # CHANGE USER
# USER $NB_USER
# ENV PYTHONPATH $PYTHONPATH:$UW2_DIR

# # setup symlink for terminal convenience 
# RUN ln -s $NB_WORK /home/$NB_USER/

# # create a volume
# VOLUME $NB_WORK/user_data
# WORKDIR $NB_WORK

# # launch notebook
# CMD ["jupyter", "notebook", "--ip='*'", "--no-browser"]








