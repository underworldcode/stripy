#! /usr/bin/env bash


THIS_FILE=`stat -f ${BASH_SOURCE[0]} || readlink -f ${BASH_SOURCE[0]}`
THIS_DIR=`dirname $THIS_FILE`
echo "Building docs in $THIS_DIR"

cd $THIS_DIR/..

pdoc --html stripy --force -o docs
mv docs/stripy/*.html docs
rmdir docs/stripy
