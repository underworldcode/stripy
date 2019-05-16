"""
Copyright 2017 Louis Moresi, Ben Mather

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


class square_border(_cartesian.Triangulation):
    """
    Square border to add to any mesh.
    Creates a set of points around the extent of the domain.
    """
    def __init__(self, extent, spacingX, spacingY, refinement_levels=0, tree=False):

        xmin, xmax, ymin, ymax = extent

        for r in range(0, refinement_levels):
            spacingX = 0.5*spacingX
            spacingY = 0.5*spacingY

        xcoords = np.arange(xmin, xmax+spacingX, spacingX)
        ycoords = np.arange(ymin, ymax+spacingY, spacingY)

        nx, ny = xcoords.size, ycoords.size

        vx = np.ones(nx-2)
        vy = np.ones(ny-2)

        x = np.concatenate([xcoords[1:-1], xcoords[1:-1], vy*xmin, vy*xmax], axis=0)
        y = np.concatenate([vx*ymin, vx*ymax, ycoords[1:-1], ycoords[1:-1]], axis=0)

        print (x.shape, y.shape)

        xy = np.column_stack([x, y])

        # Randomise the point order for triangulation efficiency
        np.random.shuffle(xy)

        # make sure the first points are not collinear
        xy0 = np.array([[xcoords[0] , ycoords[0] ], \
                        [xcoords[-1], ycoords[-1]], \
                        [xcoords[-1], ycoords[0] ], \
                        [xcoords[0] , ycoords[-1]]])

        xy = np.vstack([xy0, xy])

        super(square_border, self).__init__(xy[:,0], xy[:,1], permute=False, tree=tree)

        return


class square_mesh(_cartesian.Triangulation):
    """
    A square mesh where points are uniformly populated
    along x and y directions defined by extent.
    """
    def __init__(self, extent, spacingX, spacingY, random_scale=0.0, refinement_levels=0, tree=False):

        xmin, xmax, ymin, ymax = extent

        xcoords = np.arange(xmin, xmax+spacingX, spacingX)
        ycoords = np.arange(ymin, ymax+spacingY, spacingY)

        x, y = np.meshgrid(xcoords, ycoords)

        interior_mask = np.zeros(x.shape, dtype=bool)
        interior_mask[1:-1,1:-1] = True
        interior_mask = interior_mask.ravel()

        x = x.ravel()
        y = y.ravel()

        nsamples = interior_mask.sum()
        xscale = random_scale*spacingX
        yscale = random_scale*spacingY

        x[interior_mask] += xscale * (0.5 - np.random.rand(nsamples))
        y[interior_mask] += yscale * (0.5 - np.random.rand(nsamples))

        xy = np.column_stack([x, y])

        n2 = x.size//2

        # ensure first 3 points are not collinear
        # xy[0], xy[-1] = xy[-1].copy(), xy[0].copy()
        # xy[2], xy[n2] = xy[n2].copy(), xy[2].copy()

        ## Randomise the point order for triangulation efficiency
        np.random.shuffle(xy)

        super(square_mesh, self).__init__(xy[:,0], xy[:,1], permute=False, tree=tree)

        for r in range(0, refinement_levels):
            X, Y = self.uniformly_refine_triangulation(faces=False, trisect=False)
            self._update_triangulation(X, Y)

        return


class elliptical_mesh(_cartesian.Triangulation):
    """
    An elliptical mesh where points are successively populated at an
    increasing radius from the midpoint of the extent.
    """
    def __init__(self, extent, spacingX, spacingY, random_scale=0.0, refinement_levels=0, tree=False):

        xmin, xmax, ymin, ymax = extent

        originX = 0.5*(xmax + xmin)
        originY = 0.5*(ymax + ymin)
        radiusX = 0.5*(xmax - xmin)
        radiusY = 0.5*(ymax - ymin)
        aspect  = 0.5*(ymax - ymin)/radiusX

        spacingXY = np.hypot(spacingX, spacingY)
        radiusXY = np.hypot(radiusX, radiusY)

        radius = np.arange(spacingXY, radiusXY, spacingXY)

        x = []
        y = []
        for r in radius:
            # nsamples must be >= 3
            nsamples = max(2.0*r/spacingXY, 3.0)
            lin_samples = int(nsamples)
            ind = np.arange(0, lin_samples)

            x_ring = originX + np.cos(2.0*np.pi/lin_samples*ind)*r
            y_ring = originY + np.sin(2.0*np.pi/lin_samples*ind)*r * aspect

            x.append(x_ring)
            y.append(y_ring)

        boundary_size = len(x[-1])

        x = np.concatenate(x, axis=0)
        y = np.concatenate(y, axis=0)

        interior_mask = np.ones(x.size, dtype=bool)
        interior_mask[-boundary_size:] = False

        nsamples = interior_mask.sum()
        xscale = random_scale*spacingX
        yscale = random_scale*spacingY

        x[interior_mask] += xscale * (0.5 - np.random.rand(nsamples))
        y[interior_mask] += yscale * (0.5 - np.random.rand(nsamples))

        super(elliptical_mesh, self).__init__(x, y, permute=False, tree=tree)

        for r in range(0, refinement_levels):
            X, Y = self.uniformly_refine_triangulation(faces=True, trisect=False)
            self._update_triangulation(X, Y)

        return


class random_mesh(_cartesian.Triangulation):
    """
    A mesh of random points. Take care if you use this is parallel
    as the location of points will not be the same on all processors.


    """
    def __init__(self, extent, number_of_points=5000, tree=False):

        xmin, xmax, ymin, ymax = extent

        deltaX = xmax - xmin
        deltaY = ymax - ymin

        samples = np.sqrt(number_of_points)
        samples_x = int(samples*deltaX/deltaY)
        samples_y = int(samples*deltaY/deltaX)

        spacingX = deltaX/samples_x
        spacingY = deltaY/samples_y

        # boundary nodes
        base_mesh = square_border(extent, spacingX, spacingY)

        xy = np.random.random((number_of_points - base_mesh.npoints,2))

        # scale to extent leaving some buffer space
        x = xmin + 0.5*spacingX + xy[:,0]*(deltaX - 0.5*spacingX)
        y = ymin + 0.5*spacingY + xy[:,1]*(deltaY - 0.5*spacingY)

        x = np.concatenate([base_mesh.x, x], axis=0)
        y = np.concatenate([base_mesh.y, y], axis=0)

        super(random_mesh, self).__init__(x, y, permute=False, tree=tree)

        return



## The following are cartesian meshes warped into 3D
## such that the original x,y directions provide a natural means
## of navigating the surface.
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
