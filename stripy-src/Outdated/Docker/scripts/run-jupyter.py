#!/usr/bin/env python

from subprocess import call
import os
import time

# We want to start a server from each www directory
# where everything was built by the site-builder script

# Make sure jupyter defaults are correct (globally)

# call("jupyter nbextension enable hide_input/main", shell=True)
# call("jupyter nbextension enable rubberband/main", shell=True)
# call("jupyter nbextension enable exercise/main",   shell=True)

# This could be automated, but I am not sure how well the number of
# servers will scale ... so leave at 8 ... and hand build

# We want to start the server from the _site directory
# where everything was built by the docker-site-builder script

dir      = "Notebooks"
password = ""
port     = 8080

call( "cd {:s} && nohup jupyter notebook --port={:d} --ip='*' --no-browser \
       --NotebookApp.token={:s}  &".format(dir, port, password), shell=True )

# Don't exit

while True:
    time.sleep(10)


## --NotebookApp.default_url=/files/index.html