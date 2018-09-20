#! /usr/bin/env bash

# Template file make these substitutions:
# IMAGENAME="???"
# PROJ_NAME="???"

stripy-docker-help(){
  cat << *DOCUMENTATION*

  Installing bash functions for stripy docker shortcuts
  =====================================================

  0) Help !!

  stripy-docker-help

  1) Open the examples and serve on port 8899 by default

  stripy-docker-browse [port:-8899]

      stripy-docker-browse 8899
      python -m webbrowser -t "http://localhost:8899" # if you have a standard python

  This version starts a docker container in detached / restart mode which
  persists and will need to be manually shutdown through docker commands.
  This command requires a name for the running container which helps avoid
  multiple, clashing instances, and helps manage the processes.

  stripy-docker-serve [port:-8899]

      stripy-docker-serve name 8899
      python -m webbrowser -t "http://localhost:8899" # if you have a standard python


  2) Run a script with the docker installation version of python and all relevant
     modules/data pre-installed.

  stripy-docker-sh commands

    stripy-docker-sh python extract_values.py
                                              # Runs the docker version of python
                                              # on the local script

    stripy-docker-sh install_examples
                                              # Runs the built in command to install notebook
                                              # examples in the current local directory

  3) Run a specific notebook to access via port 8899 or browse in current
     directory if blank. Server runs on given port [or 8899 by default]

  stripy-docker-nb notebook_name.ipynb [port:-8899]

    stripy-docker-nb some-notebook.ipynb 8899
    python -m webbrowser -t "http://localhost:8899" # if you have a standard python

  4) Get into the docker image and poke around on the command line interactively

  stripy-docker-terminal

*DOCUMENTATION*
}

# Open the default version of the docker to browse the examples
stripy-docker-browse(){
  PORT=${1:-8899};
  echo "Navigate to http://localhost:$PORT to view the examples, ^\ when done";
  docker run --rm -p $PORT:8888 -v ${PWD}:/home/jovyan/external --rm lmoresi/stripy:0.7;
}

# Open the default version of the docker to serve examples (persistent)
stripy-docker-serve(){
  NAME=${1:stripy}
  PORT=${2:-8899};
  echo "Check status: docker ps | grep $NAME "
  echo "Manage:       docker stop/start/restart $NAME"
  docker run -d --restart unless-stopped -v ${PWD}:/home/jovyan/external --name $NAME -p $PORT:8888 lmoresi/stripy:0.7;
}

#
stripy-docker-sh(){
    docker run -v ${PWD}:/home/jovyan --rm lmoresi/stripy:0.7 $* ;
}

#
stripy-docker-terminal(){
    docker run -it --rm lmoresi/stripy:0.7 bash ;
}

stripy-docker-nb(){
    PORT=${2:-8899};  # default to 8899 if second argument not given
    echo "Navigate to http://localhost:$PORT to view, ^\ when done";
    docker run -v ${PWD}:/home/jovyan -p $PORT:8888 --rm lmoresi/stripy:0.7 jupyter-notebook --no-browser --ip="*" --NotebookApp.token='' --NotebookApp.open_browser=False --NotebookApp.default_url=/tree/"$1";
}

stripy-docker-help
