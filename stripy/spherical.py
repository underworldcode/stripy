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
from . import _stripack
from . import _ssrfpack
import numpy as np

try: range = xrange
except: pass

_ier_codes = {0:  "no errors were encountered.",
              -1: "N < 3 on input.",
              -2: "the first three nodes are collinear.",
              -3: "duplicate nodes were encountered.",
              -4: "an error flag was returned by a call to SWAP in ADDNOD.\n \
                   This is an internal error and should be reported to the programmer.",
              'L':"nodes L and M coincide for some M > L.\n \
                   The linked list represents a triangulation of nodes 1 to M-1 in this case.",
              1: "NCC, N, NROW, or an LCC entry is outside its valid range on input.",
              2: "the triangulation data structure (LIST,LPTR,LEND) is invalid.",
              'K': 'NPTS(K) is not a valid index in the range 1 to N.'}


class sTriangulation(object):
    """
    Define a Delaunay triangulation for given points on a sphere
    where lons and lats are 1D numpy arrays of equal length.

    Algorithm
    ---------
     R. J. Renka (1997), Algorithm 772; STRIPACK: Delaunay triangulation
     and Voronoi diagram on the surface of a sphere"
     ACM Trans. Math. Softw., Vol. 23, No. 3, pp 416-434
     doi:10.1145/275323.275329

    Parameters
    ----------
     lons : 1D array of longitudinal coordinates in radians
     lats : 1D array of latitudinal coordinates in radians

    Attributes
    ----------
     lons : array of floats, shape (n,)
        stored longitudinal coordinates on a sphere
     lats : array of floats, shape (n,)
        stored latitudinal coordinates on a sphere
     x : array of floats, shape (n,)
        stored Cartesian x coordinates from input
     y : array of floats, shape (n,)
        stored Cartesian y coordinates from input
     z : array of floats, shape (n,)
        stored Cartesian z coordinates from input
     simplices : array of ints, shape (nsimplex, 3)
        indices of the points forming the simplices in the triangulation
        points are ordered anticlockwise
     lst : array of ints, shape (6n-12,)
        nodal indices with lptr and lend, define the triangulation as a set of N
        adjacency lists; counterclockwise-ordered sequences of neighboring nodes
        such that the first and last neighbors of a boundary node are boundary
        nodes (the first neighbor of an interior node is arbitrary).  In order to
        distinguish between interior and boundary nodes, the last neighbor of
        each boundary node is represented by the negative of its index.
        The indices are 1-based (as in Fortran), not zero based (as in python).
     lptr : array of ints, shape (6n-12),)
        set of pointers in one-to-one correspondence with the elements of lst.
        lst(lptr(i)) indexes the node which follows lst(i) in cyclical
        counterclockwise order (the first neighbor follows the last neighbor).
        The indices are 1-based (as in Fortran), not zero based (as in python).
     lend : array of ints, shape (n,)
        N pointers to adjacency lists.
        lend(k) points to the last neighbor of node K.
        lst(lend(K)) < 0 if and only if K is a boundary node.
        The indices are 1-based (as in Fortran), not zero based (as in python).

    Notes
    -----
     Provided the nodes are randomly ordered, the algorithm
     has expected time complexity O(N*log(N)) for most nodal
     distributions.  Note, however, that the complexity may be
     as high as O(N**2) if, for example, the nodes are ordered
     on increasing latitude.
    """
    def __init__(self, lons=None, lats=None, refinement_levels=0):

        # lons, lats = self._check_integrity(lons, lats)

        self._update_triangulation(lons, lats)

        return


    def _update_triangulation(self, lons, lats):

        npoints = lons.size

        # compute cartesian coords on unit sphere.

        x, y, z = _stripack.trans(lats, lons)
        lst, lptr, lend, ierr = _stripack.trmesh(x, y, z)

        if ierr != 0:
            raise ValueError('ierr={} in trmesh\n{}'.format(ierr, _ier_codes[ierr]))

        self.lons = lons
        self.lats = lats
        self.npoints = npoints
        self.x = x
        self.y = y
        self.z = z
        self.points = np.column_stack([x, y, z])
        self.lst = lst
        self.lptr = lptr
        self.lend = lend

        # Convert a triangulation to a triangle list form (human readable)
        # Uses an optimised version of trlist that returns triangles
        # without neighbours or arc indices
        nt, ltri, ierr = _stripack.trlist2(lst, lptr, lend)

        if ierr != 0:
            raise ValueError('ierr={} in trlist2\n{}'.format(ierr, _ier_codes[ierr]))

        # extract triangle list and convert to zero-based ordering
        self.simplices = ltri.T[:nt] - 1

        return



    def _check_integrity(self, lons, lats):
        """
        Ensure lons and lats are:
         - 1D numpy arrays
         - equal size
         - within the appropriate range in radians
        """
        lons = np.array(lons).ravel()
        lats = np.array(lats).ravel()

        if len(lons.shape) != 1 or len(lats.shape) != 1:
            raise ValueError('lons and lats must be 1D')
        if lats.size != lons.size:
            raise ValueError('lons and lats must have same length')
        if (np.abs(lons)).max() > 2.*np.pi:
            raise ValueError("lons must be in radians (-2*pi <= lon <= 2*pi)")
        if (np.abs(lats)).max() > 0.5*np.pi:
            raise ValueError("lats must be in radians (-pi/2 <= lat <= pi/2)")
        return lons, lats


    def interpolate(self, lons, lats, zdata, order=1):
        """
        Base class to handle nearest neighbour, linear, and cubic interpolation.
        Given a triangulation of a set of nodes on the unit sphere, along with data
        values at the nodes, this method interpolates (or extrapolates) the value
        at a given longitude and latitude.

        Parameters
        ----------
         lons : float / array of floats, shape (l,)
            longitudinal coordinate(s) on the sphere
         lats : float / array of floats, shape (l,)
            latitudinal coordinate(s) on the sphere
         zdata : array of floats, shape (n,)
            value at each point in the triangulation
            must be the same size of the mesh
         order : int (default=1)
            order of the interpolatory function used
             0 = nearest-neighbour
             1 = linear
             3 = cubic

        Returns
        -------
         zi : float / array of floats, shape (l,)
            interpolates value(s) at (lons, lats)
        """
        lons, lats = self._check_integrity(lons, lats)

        if order not in [0,1,3]:
            raise ValueError("order must be 0, 1, or 3")
        if zdata.size != self.npoints:
            raise ValueError("data must be of size {}".format(self.npoints))

        zi, ierr = _stripack.interp_n(order, lats, lons,\
                                      self.x, self.y, self.z, zdata,\
                                      self.lst, self.lptr, self.lend)

        if ierr != 0:
            raise ValueError('ierr={} in interp_n\n{}'.format(ierr, _ier_codes[ierr]))

        return zi


    def interpolate_nearest(self, lons, lats, data):
        """
        Interpolate using nearest-neighbour approximation
        Returns the same as interpolate(lons,lats,data,order=0)
        """
        return self.interpolate(lons, lats, data, order=0)

    def interpolate_linear(self, lons, lats, data):
        """
        Interpolate using nearest-neighbour approximation
        Returns the same as interpolate(lons,lats,data,order=1)
        """
        return self.interpolate(lons, lats, data, order=1)

    def interpolate_cubic(self, lons, lats, data):
        """
        Interpolate using nearest-neighbour approximation
        Returns the same as interpolate(lons,lats,data,order=3)
        """
        return self.interpolate(lons, lats, data, order=3)


    def nearest_neighbour(self, lon, lat):
        """
        Get the index of the nearest neighbour to a point (lon,lat)
        and return the squared distance between (lon,lat) and
        its nearest neighbour.

        Notes
        -----
         Faster searches can be obtained using a KDTree.
         Store all x,y and z coordinates in a (c)KDTree, then query
         a set of points to find their nearest neighbours.

         The KDTree will fail if the Euclidian distance to some node is
         shorter than the great circle distance to the near neighbour
         on the surface. scipy's KDTree will also return the Euclidian
         distance whereas this routine returns the great circle distance.
        """

        # translate to unit sphere
        xi, yi, zi = _stripack.trans(lat, lon)

        # i is the node at which we start the search
        # the closest x coordinate is a good place
        i = ((self.x - xi)**2).argmin() + 1
        idx, d = _stripack.nearnd((xi,yi,zi), self.x, self.y, self.z, self.lst, self.lptr, self.lend, i)
        return idx - 1, d

    def nearest_vertex(self, lons, lats):
        """
        Locate the index of the nearest vertex to points (lons, lats)
        and return the squared great circle distance between (lons,lats) and
        each nearest neighbour.

        Notes
        -----
         Faster searches can be obtained using a KDTree.
         Store all x,y and z coordinates in a (c)KDTree, then query
         a set of points to find their nearest neighbours.

         The KDTree will fail if the Euclidian distance to some node is
         shorter than the great circle distance to the near neighbour
         on the surface. scipy's KDTree will also return the Euclidian
         distance whereas this routine returns the great circle distance.
        """

        # translate to unit sphere
        xi = np.array(_stripack.trans(lats, lons))
        idx = np.empty_like(xi[0,:], dtype=np.int)
        dist = np.empty_like(xi[0,:], dtype=np.float)

        for pt in range(0, xi.shape[1]):
            xi0 = xi[:,pt]

            # i is the node at which we start the search
            # the closest x coordinate is a good place
            i = ((self.x - xi0[0])**2).argmin() + 1

            idx0, d0 = _stripack.nearnd((xi0[0],xi0[1],xi0[2]), self.x, self.y, self.z, self.lst, self.lptr, self.lend, i)
            idx[pt]  = idx0 - 1
            dist[pt] = d0

        return idx, dist


    def gradient(self, f, nit=3, tol=1e-3, guarantee_convergence=False):
        """
        Return the gradient of an array.

        The method consists of minimizing a quadratic functional Q(G) over
        gradient vectors, where Q is an approximation to the linearized
        curvature over the triangulation of a C-1 bivariate function F(x,y)
        which interpolates the nodal values and gradients.

        Parameters
        ----------
         f : array of floats, shape (n,)
            field over which to evaluate the gradient
         nit: int (default: 3)
            number of iterations to reach a convergence tolerance, tol
            nit >= 1
         tol: float (default: 1e-3)
            maximum change in gradient between iterations.
            convergence is reached when this condition is met.

        Returns
        -------
         dfdx : array of floats, shape (n,)
            derivative of f in the x direction
         dfdy : array of floats, shape (n,)
            derivative of f in the y direction
         dfdz : array of floats, shape (n,)
            derivative of f in the z direction

        Notes
        -----
         For SIGMA = 0, optimal efficiency was achieved in testing with
         tol = 0, and nit = 3 or 4.

         The restriction of F to an arc of the triangulation is taken to be
         the Hermite interpolatory tension spline defined by the data values
         and tangential gradient components at the endpoints of the arc, and
         Q is the sum over the triangulation arcs, excluding interior
         constraint arcs, of the linearized curvatures of F along the arcs --
         the integrals over the arcs of D2F(T)**2, where D2F(T) is the second
         derivative of F with respect to distance T along the arc.
        """
        if f.size != self.npoints:
            raise ValueError('f should be the same size as mesh')

        gradient = np.zeros((3,self.npoints), order='F', dtype=np.float32)
        sigma = 0
        iflgs = 0

        ierr = 1
        while ierr == 1:
            ierr = _ssrfpack.gradg(self.x, self.y, self.z, f, self.lst, self.lptr, self.lend,\
                                   iflgs, sigma, nit, tol, gradient)
            if not guarantee_convergence:
                break

        if ierr < 0:
            raise ValueError('ierr={} in gradg\n{}'.format(ierr, _ier_codes[ierr]))

        return gradient[0], gradient[1], gradient[2]


    def transform_to_spherical(self, dfdx, dfdy, dfdz):
        """
        Transform the derivatives of f in the x,y,z directions into spherical
        derivatives.

        Arguments
        ---------
         dfdx : array of floats, shape (n,)
            derivative of f in the x direction
         dfdy : array of floats, shape (n,)
            derivative of f in the y direction
         dfdz : array of floats, shape (n,)
            derivative of f in the z direction

        Returns
        -------
         dfdlons : array of floats, shape (n,)
            derivatives of f w.r.t. longitude
         dfdlats : array of floats, shape (n,)
            derivatives of f w.r.t. latitude

        """
        cos_lons = np.cos(self.lons)
        sin_lons = np.sin(self.lons)

        cos_lats = np.cos(self.lats)
        sin_lats = np.sin(self.lats)

        dxdlons = -sin_lons*sin_lats
        dxdlats =  cos_lons*cos_lats

        dydlons =  cos_lons*sin_lats
        dydlats =  sin_lons*cos_lats

        dzdlons = -sin_lons
        dzdlats = 0.0

        # chain rule
        dfdlons = dfdx*dxdlons + dfdy*dydlons + dfdz*dzdlons
        dfdlats = dfdx*dxdlats + dfdy*dydlats

        return dfdlons, dfdlats


    def smoothing(self, f, w, sm, smtol, gstol):
        """
        Smooths a surface f by choosing nodal function values and gradients to
        minimize the linearized curvature of F subject to a bound on the
        deviation from the data values. This is more appropriate than interpolation
        when significant errors are present in the data.

        Parameters
        ----------
         f : array of floats, shape (n,)
            field to apply smoothing on
         w : array of floats, shape (n,)
            weights associated with data value in f
            w[i] = 1/sigma_f^2 is a good rule of thumb.
         sm : float
            positive parameter specifying an upper bound on Q2(f).
            generally n-sqrt(2n) <= sm <= n+sqrt(2n)
         smtol : float
            specifies relative error in satisfying the constraint
            sm(1-smtol) <= Q2 <= sm(1+smtol) between 0 and 1.
         gstol : float
            tolerance for convergence.
            gstol = 0.05*mean(sigma_f)^2 is a good rule of thumb.

        Returns
        -------
         f_smooth : array of floats, shape (n,)
            smoothed version of f
         (dfdx, dfdy, dfdz) : tuple of floats, tuple of 3 shape (n,) arrays
            first derivatives of f_smooth in the x, y, and z directions
        """
        if f.size != self.npoints or f.size != w.size:
            raise ValueError('f and w should be the same size as mesh')

        sigma = 0
        iflgs = 0
        prnt = -1

        f_smooth, df, ierr = _srfpack.smsurf(self.x, self.y, self.z, f, self.lst, self.lptr, self.lend,\
                                             iflgs, sigma, w, sm, smtol, gstol, prnt)

        if ierr < 0:
            raise ValueError('ierr={} in gradg\n{}'.format(ierr, _ier_codes[ierr]))
        if ierr == 1:
            raise RuntimeWarning("No errors were encountered but the constraint is not active --\n\
                  F, FX, and FY are the values and partials of a linear function \
                  which minimizes Q2(F), and Q1 = 0.")
        if ierr == 2:
            raise RuntimeWarning("if the constraint could not be satisfied to within SMTOL\
                  due to ill-conditioned linear systems.")

        return f_smooth, (df[0], df[1], df[2])


    def tri_area(self, lons, lats):
        """
        Calculate the area enclosed by 3 points on the unit sphere.

        Parameters
        ----------
         lons : array of floats, shape (3)
            longitudinal coordinates in radians
         lats : array of floats, shape (3)
            latitudinal coordinates in radians

        Returns
        -------
         area : float
            area of triangle on the unit sphere

        """
        lons, lats = self._check_integrity(lons, lats)

        # translate to unit sphere
        x, y, z = _stripack.trans(lats, lons)

        # compute area
        area = _stripack.areas(x, y, z)

        return area


    def areas(self):
        """
        Compute the area each triangle within the triangulation of points
        on the unit sphere.

        Returns
        -------
         area : array of floats, shape (nt,)
            area of each triangle in self.simplices where nt
            is the number of triangles.

        Notes
        -----
         This uses a Fortran 90 subroutine that wraps the AREA function
         to iterate over many points.
        """

        return _stripack.triareas(self.x, self.y, self.z, self.simplices.T+1)


##
##  Better not to have pyproj as dependency (and we assume a sphere in this module)
##
    # def sTriangulation_midpoints_pyproj(self):
    #
    #     interpolator = self
    #
    #     import pyproj
    #
    #     lst  = interpolator.lst
    #     lend = interpolator.lend
    #     lptr = interpolator.lptr
    #
    #     g = pyproj.Geod(ellps='WGS84')
    #
    #     midlon_array = np.ones((len(lptr))) * -99999.0
    #     midlat_array = np.ones((len(lptr))) * -99999.0
    #
    #     lonv1 = interpolator.lons
    #     latv1 = interpolator.lats
    #
    #     for i in range(0,len(lptr),1):
    #         n1 = lst[i]-1
    #         n2 = lst[lptr[i]-1]-1
    #         if n1 < n2:
    #             midlonlat, = g.npts(lonv1[n1],latv1[n1],lonv1[n2],latv1[n2], 1 , radians=True )
    #             midlon_array[i] = midlonlat[0]
    #             midlat_array[i] = midlonlat[1]
    #
    #     valid_points =  np.where(midlon_array != -99999.0 )
    #
    #     midlon_array = midlon_array[valid_points[0]]
    #     midlat_array = midlat_array[valid_points[0]]
    #
    #     return midlon_array, midlat_array


    def segment_midpoints(self):
        """
        Identify the midpoints of every line segment in the triangulation
        """

        lst  = self.lst
        lend = self.lend
        lptr = self.lptr

        segments_array = np.empty((len(lptr),2),dtype=np.int)
        segments_array[:,0] = lst[:] - 1
        segments_array[:,1] = lst[lptr[:]-1] - 1

        valid = np.where(segments_array[:,0] < segments_array[:,1])[0]
        segments = segments_array[valid,:]

        mids = (self.points[segments[:,0]] + self.points[segments[:,1]]) * 0.5
        mids /= np.sqrt(mids[:,0]**2 + mids[:,1]**2 + mids[:,2]**2 ).reshape(-1,1)


        midlls = self.xyz2lonlat(mids[:,0], mids[:,1], mids[:,2])

        return midlls

    def segment_tripoints(self, ratio=0.33333):
        """
        Identify the trisection points of every line segment in the triangulation
        """


        lst  = self.lst
        lend = self.lend
        lptr = self.lptr


        segments_array = np.empty((len(lptr),2),dtype=np.int)
        segments_array[:,0] = lst[:] - 1
        segments_array[:,1] = lst[lptr[:]-1] - 1

        valid = np.where(segments_array[:,0] < segments_array[:,1])[0]
        segments = segments_array[valid,:]

        mids1 = ratio * self.points[segments[:,0]] + (1.0-ratio) * self.points[segments[:,1]]
        mids1 /= np.sqrt(mids1[:,0]**2 + mids1[:,1]**2 + mids1[:,2]**2 ).reshape(-1,1)

        mids2 = (1.0-ratio) *  self.points[segments[:,0]] + ratio * self.points[segments[:,1]]
        mids2 /= np.sqrt(mids2[:,0]**2 + mids2[:,1]**2 + mids2[:,2]**2 ).reshape(-1,1)

        mids = np.vstack((mids1,mids2))

        midlls = self.xyz2lonlat(mids[:,0], mids[:,1], mids[:,2])

        return midlls


    def face_midpoints(self):
        """
        Identify the centroid of every simplex in the triangulation
        """

        mids = self.points[self.simplices].mean(axis=1)
        mids /= np.sqrt(mids[:,0]**2 + mids[:,1]**2 + mids[:,2]**2 ).reshape(-1,1)

        midlls = self.xyz2lonlat(mids[:,0], mids[:,1], mids[:,2])

        return midlls

    def lons_map_to_wrapped(self, lon):

        lons = np.array(lon)
        lons = np.mod(lon+np.pi, 2*np.pi) - np.pi

        return lons


    def lonlat2xyz(self,lon, lat):
        """
        Convert lon / lat (radians) for the spherical triangulation into x,y,z
        on the unit sphere
        """

        lons = np.array(lon)
        lats = np.array(lat)

        xs = np.cos(lats) * np.cos(lons)
        ys = np.cos(lats) * np.sin(lons)
        zs = np.sin(lats)

        return xs, ys, zs

    def xyz2lonlat(self, x,y,z):
        """
        Convert x,y,z representation of points *on the unit sphere* of the
        spherical triangulation to lon / lat (radians).

        Note - no check is made here that (x,y,z) are unit vectors
        """

        xs = np.array(x)
        ys = np.array(y)
        zs = np.array(z)

        lons = np.arctan2(y, x)
        lats = np.arcsin(z)

        return lons, lats

    def _add_spherical_midpoints(self):

        midlon_array, midlat_array = self.segment_midpoints()

        lonv2 = np.concatenate((self.lons, midlon_array), axis=0)
        latv2 = np.concatenate((self.lats, midlat_array), axis=0)

        return lonv2, latv2

    def _add_spherical_tripoints(self, ratio=0.333333):

        midlon_array, midlat_array = self.segment_tripoints(ratio=ratio)

        lonv2 = np.concatenate((self.lons, midlon_array), axis=0)
        latv2 = np.concatenate((self.lats, midlat_array), axis=0)

        return lonv2, latv2

    def _add_face_centroids(self):

        facelon_array, facelat_array = self.face_midpoints()

        lonv2 = np.concatenate((self.lons, facelon_array), axis=0)
        latv2 = np.concatenate((self.lats, facelat_array), axis=0)

        return lonv2, latv2


    def uniformly_refine_triangulation(self, faces=False, trisect=False):
        """
        return a refined triangulation obtained by bisection of all edges
        in the triangulation
        """

        if faces:
            lonv1, latv1 = self._add_face_centroids()

        else:
            if not trisect:
                lonv1, latv1 = self._add_spherical_midpoints()
            else:
                lonv1, latv1 = self._add_spherical_tripoints(ratio=0.333333)


        return lonv1, latv1

    def join(self, t2, unique=False):
        """
        Join this triangulation with another. If the points are known to have no duplicates, then
        set unique=False to skip the testing and duplicate removal
        """

        lonv1 = np.concatenate((self.lons, t2.lons), axis=0)
        latv1 = np.concatenate((self.lats, t2.lats), axis=0)

        ## remove any duplicates

        if not unique:
            a = np.ascontiguousarray(np.vstack((lonv1, latv1)).T)
            unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
            llunique = unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

            lonv1 = llunique[:,0]
            latv1 = llunique[:,1]

        return lonv1, latv1
