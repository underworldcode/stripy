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


class square_mesh(_cartesian.Triangulation):
    """
    A square mesh where points are uniformly populated
    along x and y directions defined by extent
    """
    def __init__(self, extent, spacingX, spacingY, number_of_points, refinement_levels=0):

        xmin, xmax, ymin, ymax = extent

        nsamples = np.sqrt(number_of_points)
        lin_samples = int(nsamples)

        tiX = np.linspace(xmin + 0.75 * spacingX, xmax - 0.75 * spacingX, lin_samples)
        tiY = np.linspace(ymin + 0.75 * spacingY, ymax - 0.75 * spacingY, lin_samples)

        x, y = np.meshgrid(tiX, tiY)

        x = x.ravel()
        y = y.ravel()

        xscale = (x.max()-x.min()) / (2.0 * nsamples)
        yscale = (y.max()-y.min()) / (2.0 * nsamples)

        x += xscale * (0.5 - np.random.rand(x.size))
        y += yscale * (0.5 - np.random.rand(y.size))

        super(square_mesh, self).__init__(x, y)

        for r in range(0, refinement_levels):
            x, y = self.uniformly_refine_triangulation(faces=False, trisect=False)
            self._update_triangulation(x, y)

        return


class elliptical_mesh(_cartesian.Triangulation):
    """
    An elliptical mesh where points are populated at a radius
    of the midpoint defined by extent
    """
    def __init__(self, extent, spacingX, spacingY, number_of_points, refinement_levels=0):

        xmin, xmax, ymin, ymax = extent

        originX = 0.5 * (xmax + xmin)
        originY = 0.5 * (ymax + ymin)
        radiusX = 0.5 * (xmax - xmin)
        aspect  = 0.5 * (ymax - ymin) / radiusX

        nsamples = np.sqrt(number_of_points)
        lin_samples = int(nsamples)

        tiX = np.linspace(xmin + 0.75 * spacingX, xmax - 0.75 * spacingX, lin_samples)
        tiY = np.linspace(ymin + 0.75 * spacingY, ymax - 0.75 * spacingY, lin_samples)

        x, y = np.meshgrid(tiX, tiY)

        x = np.reshape(x,len(x)*len(x[0]))
        y = np.reshape(y,len(y)*len(y[0]))

        xscale = (x.max()-x.min()) / (2.0 * nsamples)
        yscale = (y.max()-y.min()) / (2.0 * nsamples)

        x += xscale * (0.5 - np.random.rand(len(x)))
        y += yscale * (0.5 - np.random.rand(len(y)))

        mask = np.where( (x**2 + y**2 / aspect**2) < (radiusX-0.5*spacingX)**2 )

        X = x[mask]
        Y = y[mask]

        super(elliptical_mesh, self).__init__(X, Y)

        for r in range(0, refinement_levels):
            X, Y = self.uniformly_refine_triangulation(faces=False, trisect=False)
            self._update_triangulation(X, Y)

        return


class random_mesh(_cartesian.Triangulation):
    """
    A mesh of random points. Take care if you use this is parallel
    as the location of points will not be the same on all processors
    """
    def __init__(self, extent, number_of_points=5000):

        xmin, xmax, ymin, ymax = extent

        x = np.random.uniform(xmin, xmax, size=number_of_points)
        y = np.random.uniform(ymin, ymax, size=number_of_points)

        super(random_mesh, self).__init__(x, y)

        return