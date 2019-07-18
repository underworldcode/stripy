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
from . import spherical as _spherical
import numpy as np




class icosahedral_mesh(_spherical.sTriangulation):
    """
    An icosahedral triangulated mesh based on the sTriangulation class.
    """

    def __init__(self, refinement_levels=0, include_face_points=False, trisection=False, tree=False):

        mid_lat = np.degrees(np.arctan(0.5)) # 26.56505117707799 degrees
        vertices_LatLonDeg = np.array(
                                       [[ 90,         0.0 ],
                                        [ mid_lat,    0.0 ],
                                        [-mid_lat,   36.0 ],
                                        [ mid_lat,   72.0 ],
                                        [-mid_lat,  108.0 ],
                                        [ mid_lat,  144.0 ],
                                        [-mid_lat,  180.0 ],
                                        [ mid_lat,  360.0-72.0 ],
                                        [-mid_lat,  360.0-36.0 ],
                                        [ mid_lat,  360.0-144.0 ],
                                        [-mid_lat,  360.0-108.0 ],
                                        [-90,         0.0 ]])


        vertices_lat = np.radians(vertices_LatLonDeg.T[0])
        vertices_lon = np.radians(vertices_LatLonDeg.T[1])

        super(icosahedral_mesh, self).__init__(lons=vertices_lon, lats=vertices_lat, permute=False, tree=tree)

        if include_face_points:
            lons, lats = self.uniformly_refine_triangulation(faces=True)
            self._update_triangulation(lons,lats)

        for r in range(0,refinement_levels):
            lons, lats = self.uniformly_refine_triangulation(faces=False, trisect=trisection)
            self._update_triangulation(lons,lats)

        return


class octahedral_mesh(_spherical.sTriangulation):
    """
    An octahedral triangulated mesh based on the sTriangulation class
    """

    def __init__(self, refinement_levels=0, include_face_points=False, tree=False):

        vertices_LatLonDeg = np.array(
                                       [[  0.00,    0.0],
                                        [ 90.0,     0.0],
                                        [  0.00,   90.0],
                                        [  0.00,  180.0],
                                        [ -0.00,  270.0],
                                        [-90.0,    180.0]]
                                    )


        vertices_lat = np.radians(vertices_LatLonDeg.T[0])
        vertices_lon = np.radians(vertices_LatLonDeg.T[1])

        super(octahedral_mesh, self).__init__(lons=vertices_lon, lats=vertices_lat, permute=False, tree=tree)

        if include_face_points:
            lons, lats = self.uniformly_refine_triangulation(faces=True)
            self._update_triangulation(lons,lats)


        for r in range(0,refinement_levels):
            lons, lats = self.uniformly_refine_triangulation(faces=False, trisect=False )
            self._update_triangulation(lons,lats)

        return


class triangulated_cube_mesh(_spherical.sTriangulation):
    """
    An cube-based triangulated mesh based on the sTriangulation class
    """

    def __init__(self, refinement_levels=0, tree=False):

        mid_lat = np.degrees(np.arctan(np.sqrt(0.5))) # 35.264389682754654 degrees
        vertices_LatLonDeg = np.array(
                                         [[  mid_lat,   0.0],
                                          [ -mid_lat,   0.0],
                                          [  mid_lat,  90.0],
                                          [ -mid_lat,  90.0],
                                          [  mid_lat, 180.0],
                                          [ -mid_lat, 180.0],
                                          [  mid_lat, -90.0],
                                          [ -mid_lat, -90.0]   ]
                                    )


        vertices_lat = np.radians(vertices_LatLonDeg.T[0])
        vertices_lon = np.radians(vertices_LatLonDeg.T[1])

        super(triangulated_cube_mesh, self).__init__(lons=vertices_lon, lats=vertices_lat, permute=False, tree=tree)

        for r in range(0,refinement_levels):
            lons, lats = self.uniformly_refine_triangulation(faces=False, trisect=False)
            self._update_triangulation(lons,lats)

        return


class triangulated_soccerball_mesh(_spherical.sTriangulation):
    """
    This mesh is inspired by the C60 molecule and the soccerball - a truncated
    icosahedron with mid points added to all pentagon and hexagon faces to create
    a uniform triangulation.
    """

    def __init__(self, refinement_levels=0, tree=False):

        base_mesh = icosahedral_mesh(refinement_levels=1, trisection=True)
        face_mesh = icosahedral_mesh(refinement_levels=0, include_face_points=True)

        lons, lats = base_mesh.join(face_mesh)

        ll = np.vstack((lons, lats)).T

        ## Now randomise the point order
        np.random.shuffle(ll)

        super(triangulated_soccerball_mesh, self).__init__(lons=ll[:,0], lats=ll[:,1], permute=True, tree=tree)

        for r in range(0,refinement_levels):
            lons, lats = self.uniformly_refine_triangulation(faces=False, trisect=False)
            self._update_triangulation(lons,lats)

        return

class random_mesh(_spherical.sTriangulation):
    """
    A mesh of random points. Take care if you use this is parallel
    as the location of points will not be the same on all processors
    """

    def __init__(self, number_of_points=5000, tree=False):

        xyz =  np.random.random((number_of_points,3)) * 2.0 - 1.0
        xyz /= np.sqrt(xyz[:,0]**2 + xyz[:,1]**2 + xyz[:,2]**2).reshape(-1,1)

        lon,lat = _spherical.xyz2lonlat(xyz[:,0], xyz[:,1], xyz[:,2])

        super(random_mesh, self).__init__(lons=lon, lats=lat, permute=True, tree=tree)

        return


class uniform_ring_mesh(_spherical.sTriangulation):
    """
    A mesh of made of rings to create a roughly gridded, even spacing on
    the sphere. There is a small random component to prevent points lying along the
    prime meridian so this mesh should be used with caution in parallel
    """


    def __init__(self, resolution=9, refinement_levels=0, tree=False):

        # mesh of uniform-spaced rings

        eq_points = resolution
        spacing = np.pi / eq_points

        lat_values = np.linspace(-np.pi/2.0, np.pi/2.0, num=eq_points+1)[1:-1]

        lats_list = [np.array([-np.pi*0.5, np.pi*0.5])]
        lons_list = [np.array([ 0.0,       0.0      ])]

        for lat in lat_values:
            length = np.fabs(np.cos(lat)) * 2 * np.pi
            samples = int(length // spacing)
            offset = (2.0 * lat / samples)

            lons = offset + np.linspace(0.0, 2.0*np.pi, samples, endpoint=False)
            lats = np.ones_like(lons)*lat

            lons_list.append(lons)
            lats_list.append(lats)

        lons_array = np.concatenate(lons_list)
        lats_array = np.concatenate(lats_list)

        ll = np.vstack((lons_array, lats_array)).T


        ## Now randomise the point order
        np.random.shuffle(ll)

        super(uniform_ring_mesh, self).__init__(lons=ll[:,0], lats=ll[:,1], permute=False, tree=tree)

        for r in range(0,refinement_levels):
            lons, lats = self.uniformly_refine_triangulation(faces=False, trisect=False)
            self._update_triangulation(lons,lats)


        return
