#! /usr/bin/env bash

# Template file make these substitutions:
# IMAGENAME="???"
# PROJ_NAME="???"

$PROJ_NAME-docker-help(){
  cat << *DOCUMENTATION*

  Installing bash functions for $PROJ_NAME docker shortcuts
  =====================================================

  0) Help !!

  $PROJ_NAME-docker-help

  1) Open the examples and serve on port 8899 by default

  $PROJ_NAME-docker-browse [port:-8899]

      $PROJ_NAME-docker-browse 8899
      python -m webbrowser -t "http://localhost:8899" # if you have a standard python

  This version starts a docker container in detached / restart mode which
  persists and will need to be manually shutdown through docker commands.
  This command requires a name for the running container which helps avoid
  multiple, clashing instances, and helps manage the processes.

  $PROJ_NAME-docker-serve [port:-8899]

      $PROJ_NAME-docker-serve name 8899
      python -m webbrowser -t "http://localhost:8899" # if you have a standard python


  2) Run a script with the docker installation version of python and all relevant
     modules/data pre-installed.

  $PROJ_NAME-docker-sh commands

    $PROJ_NAME-docker-sh python extract_values.py
                                              # Runs the docker version of python
                                              # on the local script

    $PROJ_NAME-docker-sh install_examples
                                              # Runs the built in command to install notebook
                                              # examples in the current local directory

  3) Run a specific notebook to access via port 8899 or browse in current
     directory if blank. Server runs on given port [or 8899 by default]

  $PROJ_NAME-docker-nb notebook_name.ipynb [port:-8899]

    $PROJ_NAME-docker-nb some-notebook.ipynb 8899
    python -m webbrowser -t "http://localhost:8899" # if you have a standard python

  4) Get into the docker image and poke around on the command line interactively

  $PROJ_NAME-docker-terminal

*DOCUMENTATION*
}

# Open the default version of the docker to browse the examples
$PROJ_NAME-docker-browse(){
  PORT=${1:-8899} ;
  echo "Navigate to http://localhost:${PORT} to view the examples, ^\ when done" ;
  docker run -it --rm -p ${PORT}:8888 -v ${PWD}:/home/jovyan/external $IMAGENAME ;
}

# Open the default version of the docker to serve examples (persistent)
$PROJ_NAME-docker-serve(){
  NAME=${1:$PROJ_NAME} ;
  PORT=${2:-8899} ;
  echo "Check status: docker ps | grep $NAME " ;
  echo "Manage:       docker stop/start/restart $NAME" ;
  docker run -d --restart unless-stopped -v ${PWD}:/home/jovyan/external --name $NAME -p ${PORT}:8888 $IMAGENAME
}

#
$PROJ_NAME-docker-sh(){
    docker run -v ${PWD}:/home/jovyan --rm $IMAGENAME $* ;
}

#
$PROJ_NAME-docker-terminal(){
    docker run -it --rm $IMAGENAME bash ;
}

$PROJ_NAME-docker-nb(){
    PORT=${2:-8899};  # default to 8899 if second argument not given
    echo "Navigate to http://localhost:$PORT to view, ^\ when done";
    docker run -v ${PWD}:/home/jovyan -p $PORT:8888 --rm $IMAGENAME jupyter-notebook --no-browser --ip="*" --NotebookApp.token='' --NotebookApp.open_browser=False --NotebookApp.default_url=/tree/"$1";
}

$PROJ_NAME-docker-help
