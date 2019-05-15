"""
Copyright 2017 Louis Moresi, Ben Mather

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

from .spherical import sTriangulation
from .cartesian import Triangulation
from . import spherical_meshes
from . import cartesian_meshes
from . import documentation

## The following functions are general across sTriangulations and Triangulations

def weighted_average_to_nodes(x1, x2, data, interpolator ):
    """ Weighted average of scattered data to the nodal points
    of a triangulation using the barycentric coordinates as
    weightings.

    Parameters
    ----------
     x1, x2 : 1D arrays arrays of x,y or lon, lat (radians)
     data :   1D array of data to be lumped to the node locations
     interpolator : a stripy.Triangulation or stripy.sTriangulation object
     which defines the node locations and their triangulation

    Returns
    -------
     grid  : 1D array containing the results of the weighted average
     norm  : 1D array of the normalisation used to compute `grid`
     count : 1D int array of number of points that contribute anything to a given node

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
