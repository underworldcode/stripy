#! /usr/bin/env python

"""
Copyright 2016-2017 Louis Moresi, Ben Mather, Romain Beucher

This file is part of Quagmire.

Quagmire is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or any later version.

Quagmire is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Quagmire.  If not, see <http://www.gnu.org/licenses/>.
"""

from distutils import dir_util as _dir_util
import os

# # Update=True only overwrites files that are older (caution if updating files in the container)
# Notebooks_Path = "/home/user_data/Data"
# Destination_Path = "/home/jovyan/PYE_2019/Data"
#
# ct = _dir_util.copy_tree(Notebooks_Path, Destination_Path, preserve_mode=1, preserve_times=1, preserve_symlinks=1, update=True, verbose=True, dry_run=0)

# Update=True only overwrites files that are older (caution if updating files in the container)
Notebooks_Path = "/home/user_data/Notebooks"
Destination_Path = "/home/jovyan/STRIPY/Notebooks"

ct = _dir_util.copy_tree(Notebooks_Path, Destination_Path, preserve_mode=1, preserve_times=1, preserve_symlinks=1, update=True, verbose=True, dry_run=0)
