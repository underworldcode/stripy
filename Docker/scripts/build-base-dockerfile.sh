#!/bin/sh

# Docker build command for the base python image
# Currently we are at version 1.05

# 1.05 includes the jupyter notebook extensions but not litho1pt0 and stripy
# The latter are built at the next level up (using the modules in this repo)


set -e
cd $(dirname "$0")/../..

echo `pwd`

docker build -f Docker/Dockerfile -t lmoresi/stripy:1.05 .
