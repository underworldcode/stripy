#! /usr/bin/env bash

jupyter-book build . 

# This is best done by hand so it updates the slides even 
# if there is no work to be done in rebuilding the book 

mkdir -p _build/html/images/k3d
cp -r images/k3d _build/html/images
