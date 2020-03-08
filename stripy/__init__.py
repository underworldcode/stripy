"""
Copyright 2017-2019 Louis Moresi, Ben Mather

Stripy is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or any later version.

Stripy is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Stripy.  If not, see <http://www.gnu.org/licenses/>.

Stripy source code is available from <https://github.com/underworldcode/stripy>

"""


import os as _os
from platform import system as _system

# add '.dll' files if we are on Windows
if _system() == "Windows":
    extra_dll_dir = _os.path.join(_os.path.dirname(__file__), 'extra-dll')
    if _os.path.isdir(extra_dll_dir):
        _os.environ["PATH"] += _os.pathsep + extra_dll_dir
    extra_dll_dir = _os.path.join(_os.path.dirname(__file__), '.libs')
    if _os.path.isdir(extra_dll_dir):
        _os.environ["PATH"] += _os.pathsep + extra_dll_dir


from .spherical import sTriangulation
from .cartesian import Triangulation
from . import spherical_meshes
from . import cartesian_meshes
from . import hybrid_st_meshes
from . import documentation

## The following functions are general across sTriangulations and Triangulations

def weighted_average_to_nodes(x1, x2, data, interpolator ):
    """
    Weighted average of scattered data to the nodal points
    of a triangulation using the barycentric coordinates as
    weightings.

    Args:
        x1 : 1D array
            x,y or lon, lat (radians)
        x2 : 1D array
            x,y or lon, lat (radians)
        data : 1D array
            1D array of data to be lumped to the node locations
        interpolator : object
            a `stripy.Triangulation` or `stripy.sTriangulation` object
            which defines the node locations and their triangulation

    Returns:
        grid  : 1D array
            contains the results of the weighted average
        norm  : 1D array
            normalisation used to compute `grid`
        count : 1D int array
            number of points that contribute anything to a given node
    """

    import numpy as np

    gridded_data = np.zeros(interpolator.npoints)
    norm         = np.zeros(interpolator.npoints)
    count        = np.zeros(interpolator.npoints, dtype=np.int)

    bcc, nodes = interpolator.containing_simplex_and_bcc(x1, x2)

    # Beware vectorising the reduction operation !!

    for i in range(0, len(data)):

        grid[nodes[i][0]] += bcc[i][0] * data[i]
        grid[nodes[i][1]] += bcc[i][1] * data[i]
        grid[nodes[i][2]] += bcc[i][2] * data[i]

        norm[nodes[i][0]] += bcc[i][0]
        norm[nodes[i][1]] += bcc[i][1]
        norm[nodes[i][2]] += bcc[i][2]

        count[nodes[i][0]] += 1
        count[nodes[i][1]] += 1
        count[nodes[i][2]] += 1


    grid[np.where(norm > 0.0)] /= norm[np.where(norm > 0.0)]

    return grid, norm, count


def remove_duplicate_points(vector_tuple):
    """
    Remove duplicates rows from N equally-sized arrays
    """
    array = np.column_stack(vector_tuple)
    a = np.ascontiguousarray(array)
    unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
    b = unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))
    return list(b.T)
