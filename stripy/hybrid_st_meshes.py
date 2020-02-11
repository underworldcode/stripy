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

#!/usr/bin/python
# -*- coding: utf-8 -*-
from . import cartesian as _cartesian
import numpy as np

## The following are cartesian meshes warped into 3D
## such that the original x,y directions provide a natural metric
## for navigating the surface.
##
## These are useful for wrapping textures on surfaces
##
## Sphere (with x/y -> lon/lat) and singular at the poles
## which cannot be meshed with sTriangulation
##
## Currently: XX,YY,ZZ are the R3 coords of the points
##            SS, TT are the normalised coordinates on the manifold
##            x,y are the original point locations
##            areas etc all refer to the undistorted mesh ( ... )
##
## Also possible: cylinder / torus (less useful for visualisation, I suppose ... LM)


class warped_xy_mesh_sphere(_cartesian.Triangulation):
    """
    A lon/lat mesh on the sphere which, due to the poles, cannot be
    meshed directly with sTriangulation routines.

    This is a hybrid mesh which is calculated in a flat projection
    that corresponds to cartopy.crs.PlateCarree() but has associated
    (x,y,z) coordinates.

    The primary use case is visualisation with texture maps where the
    (s,t) coordinate system that is required on the surface has to map
    correctly for an (x,y) array of pixels.
    """

    def __init__(self, res_lon, res_lat, epsilon=0.001):

        import numpy as np

        lon = np.linspace(0.0, 2.0*np.pi, res_lon, endpoint=True)
        lat = np.linspace(epsilon,  np.pi-epsilon,  res_lat, endpoint=True)
        lons, lats = np.meshgrid(lon, lat)
        lons = lons.reshape(-1)
        lats = lats.reshape(-1)

        XX = np.cos(lons) * np.sin(lats)
        YY = np.sin(lons) * np.sin(lats)
        ZZ = np.cos(lats)

        SS = lons / (2.0*np.pi)
        TT = lats / np.pi

        super(warped_xy_mesh_sphere, self).__init__(x=SS, y=TT, permute=True, tree=False)

        self.XX = XX
        self.YY = YY
        self.ZZ = ZZ

        self.SS = SS
        self.TT = TT

        ## how should we re-define things like self.areas, self.tree ??
