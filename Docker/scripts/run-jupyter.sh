#!/usr/bin/env sh

PASS=${NB_PASSWD:-""}
OPEN=${START_NB:-""}
PORT=${NB_PORT:-8888}
#
cd $NB_DIR

jupyter-notebook --port=$NB_PORT --ip='0.0.0.0' --no-browser --allow-root \
       --NotebookApp.token=$NB_PASSWD --NotebookApp.default_url="/tree/$START_NB"

# Don't exit

while true; do
  sleep 600
done
