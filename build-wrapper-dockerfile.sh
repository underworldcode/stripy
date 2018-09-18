#!/bin/sh

# Don't forget to increment the version number if you want to keep the old stuff

IMAGENAME="lmoresi/uom-py-lavavu-notebook-bundle:1.0.7"
FROM_IMG="lmoresi/unimelb-debian-python:1.0.7"

docker build -t $IMAGENAME \
             -f Docker/WrapperDockerfile \
             --build-arg FROMIMG_ARG=$FROM_IMG \
             --build-arg IMAGENAME_ARG=$IMAGENAME \
             $PWD
