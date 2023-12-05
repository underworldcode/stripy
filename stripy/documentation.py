"""
Copyright 2017-2019 Louis Moresi, Ben Mather

This file is part of Stripy.

Stripy is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or any later version.

Stripy is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Stripy.  If not, see <http://www.gnu.org/licenses/>.
"""

import os as _os
import shutil as _shutil
from importlib import resources as _importlib_resources
from platform import python_version as _python_version


def _copytree(src, dst, symlinks=False, ignore=None):
    for item in _os.listdir(src):
        s = _os.path.join(src, item)
        d = _os.path.join(dst, item)
        if _os.path.isdir(s):
            _shutil.copytree(s, d, symlinks, ignore)
        else:
            _shutil.copy2(s, d)


def install_documentation(path="./Stripy-Notebooks"):
    """Install the example notebooks for stripy in the given location

    WARNING: If the path exists, the Notebook files will be written into the path
    and will overwrite any existing files with which they collide. The default
    path ("./Stripy-Notebooks") is chosen to make collision less likely / problematic

    The documentation for stripy is in the form of jupyter notebooks.

    Some dependencies exist for the notebooks to be useful:

       - matplotlib: for some of the diagrams
       - lavavu: interactive viewer for an alternative means of viewing cartesian_meshes
       - cartopy: for plotting map examples

    Stripy dependencies are explicitly imported into the notebooks including:

       - numpy
       - scipy (for k-d tree point location)

    """

    ## Question - overwrite or not ? shutils fails if directory exists.
    version = _python_version().split('.')
    python_more_than_v39 = int(version[0]) >=3 and int(version[1]) >= 9

    if python_more_than_v39:
        ref = _importlib_resources.files('stripy') / 'Notebooks'
        with _importlib_resources.as_file(ref) as path:
            _copytree(Notebooks_Path, path)

    else:
        import pkg_resources as _pkg_resources
        Notebooks_Path = _pkg_resources.resource_filename('stripy', 'Notebooks')
        _copytree(Notebooks_Path, path)
        # ct = _dir_util.copy_tree(Notebooks_Path, path, preserve_mode=1, preserve_times=1, preserve_symlinks=1, update=0, verbose=1, dry_run=0)

    return
