#!/usr/bin/env bash

PASS=${NB_PASSWD:-""}
OPEN=${START_NB:-""}
PORT=${NB_PORT:-8888}

cd /home/jovyan/Notebooks

jupyter-notebook --port=$NB_PORT --ip='*' --no-browser \
       --NotebookApp.token=$NB_PASSWD --NotebookApp.default_url=/tree/$OPEN

# Don't exit

while true; do
  sleep 600
done
