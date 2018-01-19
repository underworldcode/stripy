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

    def __init__(self, lons=None, lats=None, refinement_levels=0, tree=False):

        # lons, lats = self._check_integrity(lons, lats)

        self.tree = tree

        self._update_triangulation(lons, lats)

        for r in range(0,refinement_levels):
            lons, lats = self.uniformly_refine_triangulation(faces=False, trisect=False)
            self._update_triangulation(lons,lats)



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
        ## np.ndarray.sort(self.simplices, axis=1)

        ## If scipy is installed, build a KDtree to find neighbour points

        if self.tree:
            self._build_cKDtree()

        return

    def gradient_lonlat(self, data, nit=3, tol=1.0e-3, guarantee_convergence=False):
        """
        Return the lon / lat components of the gradient
        of a scalar field on the surface of the sphere.


        The method consists of minimizing a quadratic functional Q(G) over
        gradient vectors, where Q is an approximation to the linearized
        curvature over the triangulation of a C-1 bivariate function F(x,y)
        which interpolates the nodal values and gradients.

        Parameters
        ----------
         data : array of floats, shape (n,)
            field over which to evaluate the gradient
         nit: int (default: 3)
            number of iterations to reach a convergence tolerance, tol
            nit >= 1
         tol: float (default: 1e-3)
            maximum change in gradient between iterations.
            convergence is reached when this condition is met.

        Returns
        -------
         dfdlon : array of floats, shape (n,)
            derivative of f in the longitudinal direction
         dfdlat : array of floats, shape (n,)
            derivative of f in the lattitudinal direction

        Notes
        -----

        The gradient is computed via the Cartesian components using
        spherical.sTriangulation.gradient_xyz and the iteration parameters
        controling the spline interpolation are passed directly to this
        routine (See notes for gradient_xyz for more details).

        The gradient operator in this geometry is not well defined at the poles
        even if the scalar field is smooth and the Cartesian gradient is well defined.

        The routine spherical.dxyz2dlonlat is available to convert the Cartesian
        to lon/lat coordinates at any point on the unit sphere. This is helpful
        to avoid recalculation if you need both forms.
        """

        dfxs, dfys, dfzs = self.gradient_xyz(data, nit=nit, tol=tol, guarantee_convergence=guarantee_convergence)

        lons = self.lons
        lats = self.lats
        z = self.z

        dlon = -dfxs * np.cos(lats) * np.sin(lons) + dfys * np.cos(lats) * np.cos(lons) # no z dependence
        dlat = -dfxs * np.sin(lats) * np.cos(lons) - dfys * np.sin(lats) * np.sin(lons) + dfzs * np.cos(lats)

        corr = np.sqrt((1.0-z**2))
        valid = ~np.isclose(corr,0.0)

        dlon[valid] = dlon[valid] / corr[valid]

        return dlon , dlat


    def gradient_xyz(self, f, nit=3, tol=1e-3, guarantee_convergence=False):
        """
        Return the cartesian components of the gradient
        of a scalar field on the surface of the sphere.

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


##
## This one is broken - derivatives have incorrect terms and implementation of the
## gradient operator is not complete for the surface of a sphere.
##

    # def transform_to_spherical(self, dfdx, dfdy, dfdz):
    #     """
    #     Transform the Cartesian derivatives of f in the x,y,z directions into spherical
    #     derivatives.
    #
    #     Arguments
    #     ---------
    #      dfdx : array of floats, shape (n,)
    #         derivative of f in the x direction
    #      dfdy : array of floats, shape (n,)
    #         derivative of f in the y direction
    #      dfdz : array of floats, shape (n,)
    #         derivative of f in the z direction
    #
    #     Returns
    #     -------
    #      dfdlons : array of floats, shape (n,)
    #         derivatives of f w.r.t. longitude
    #      dfdlats : array of floats, shape (n,)
    #         derivatives of f w.r.t. latitude
    #
    #     """
    #     cos_lons = np.cos(self.lons)
    #     sin_lons = np.sin(self.lons)
    #
    #     cos_lats = np.cos(self.lats)
    #     sin_lats = np.sin(self.lats)
    #
    #     dxdlons = -sin_lons*sin_lats
    #     dxdlats =  cos_lons*cos_lats
    #
    #     dydlons =  cos_lons*sin_lats
    #     dydlats =  sin_lons*cos_lats
    #
    #     dzdlons = -sin_lons  # this seems wrong to me z and lat should be connected ?
    #     dzdlats = 0.0
    #
    #     # chain rule
    #     dfdlons = dfdx*dxdlons + dfdy*dydlons + dfdz*dzdlons
    #     dfdlats = dfdx*dxdlats + dfdy*dydlats
    #
    #     return dfdlons, dfdlats


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

        f_smooth, df, ierr = _ssrfpack.smsurf(self.x, self.y, self.z, f, self.lst, self.lptr, self.lend,\
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

        import warnings

        shape = lons.shape

        lons, lats = self._check_integrity(lons, lats)

        if order not in [0,1,3]:
            raise ValueError("order must be 0, 1, or 3")
        if zdata.size != self.npoints:
            raise ValueError("data must be of size {}".format(self.npoints))

        zi, zierr, ierr = _stripack.interp_n(order, lats, lons,\
                                      self.x, self.y, self.z, zdata,\
                                      self.lst, self.lptr, self.lend)

        if ierr != 0:
            print 'Warning some points may have errors - check error array\n'.format(ierr)
            zi[zierr < 0] = float('NaN')

        return zi.reshape(shape), zierr.reshape(shape)


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
        Get the index of the nearest vertex to a given point (lon,lat)
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

        Parameters
        ----------
         lons : float / array of floats, shape (l,)
            longitudinal coordinate(s) on the sphere
         lats : float / array of floats, shape (l,)
            latitudinal coordinate(s) on the sphere

        Returns
        -------

         index: array of ints - the nearest vertex to each of the supplied points

         distance: array of floats
            great circle distance (angle) on the unit sphere to the closest
            vertex identified.


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

    def containing_triangle(self, lons, lats):
        """
        Returns indices of the triangles containing lons / lats.

        Parameters
        ----------
         lons : float / array of floats, shape (l,)
            longitudinal coordinate(s) on the sphere
         lats : float / array of floats, shape (l,)
            latitudinal coordinate(s) on the sphere

        Returns
        -------

         tri_indices: array of ints, shape (l,)


        Notes
        -----
          The simplices are found as spherical.sTriangulation.simplices[tri_indices]

        """

        pts = np.array(lonlat2xyz(lons,lats)).T

        sorted_simplices = np.sort(self.simplices, axis=1)

        triangles = []
        for p in pts:
            t = _stripack.trfind(3, p, self.x, self.y, self.z, self.lst, self.lptr, self.lend )
            tri = np.sort(t[3:6])-1

            triangles.append(np.where(np.all(sorted_simplices==tri, axis=1))[0])

        return np.array(triangles).reshape(-1)

    def containing_simplex_and_bcc(self, lons, lats):
        """
        Returns the simplices containing lons / lats
        and the local barycentric, normalised coordinates.

        Notes:
        ------

        That the ordering
        of the vertices may differ from that stored in the self.simplices
        array but will still be a loop around the simplex.
        """

        pts = np.array(lonlat2xyz(lons,lats)).T

        simplices = []
        bccs = []

        for p in pts:
            t = _stripack.trfind(3, p, self.x, self.y, self.z, self.lst, self.lptr, self.lend )
            tri = t[3:6]
            bcc = t[0:3]

            simplices.append(tri)
            bccs.append(bcc)

        bcc_array = np.array(bccs)
        bcc_array /= bcc_array.sum(axis=1).reshape(-1,1)

        simplices_array = np.array(simplices)-1

        return bcc_array, simplices_array


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

    def identify_vertex_neighbours(self, vertex):
        """
        Find the neighbour-vertices in the triangulation for the given vertex
        (from the data structures of the triangulation)
        """

        lpl = self.lend[vertex-1]
        lp = lpl

        neighbours = []

        while True:
            lp = self.lptr[lp-1]
            neighbours.append(self.lst[lp-1]-1)
            if (lp == lpl):
                break

        return neighbours


    def identify_vertex_triangles(self, vertices):
        """
        Find all triangles which own any of the vertices in the list provided
        """

        triangles = []

        for vertex in np.array(vertices).reshape(-1):
            triangles.append(np.where(self.simplices == vertex)[0])

        return  np.unique(np.concatenate(triangles))



    def identify_segments(self):
        """
        Find all the segments in the triangulation and return an
        array of vertices (n1,n2) where n1 < n2
        """

        lst  = self.lst
        lend = self.lend
        lptr = self.lptr

        segments_array = np.empty((len(lptr),2),dtype=np.int)
        segments_array[:,0] = lst[:] - 1
        segments_array[:,1] = lst[lptr[:]-1] - 1

        valid = np.where(segments_array[:,0] < segments_array[:,1])[0]
        segments = segments_array[valid,:]

        return segments

    def segment_midpoints_by_vertices(self, vertices):
        """
        Add midpoints to any segment connected to the vertices in the
        list / array provided.
        """

        segments = []

        for vertex in vertices:
            neighbours = self.identify_vertex_neighbours(vertex)

            for n1 in neighbours:
                segments.append( min( tuple((vertex, n1)), tuple((n1, vertex))) )

        segs = np.array(list(set(segments)))

        new_midpoint_lonlats = self.segment_midpoints(segments=segs)

        return new_midpoint_lonlats


    def face_midpoints(self, simplices=None):
        """
        Identify the centroid of every simplex in the triangulation. If an array of
        simplices is given then the centroids of only those simplices is returned.
        """

        if type(simplices) == type(None):
            simplices = self.simplices

        mids = self.points[simplices].mean(axis=1)
        mids /= np.sqrt(mids[:,0]**2 + mids[:,1]**2 + mids[:,2]**2 ).reshape(-1,1)

        midlons, midlats = xyz2lonlat(mids[:,0], mids[:,1], mids[:,2])

        return midlons, midlats


    def segment_midpoints(self, segments=None):
        """
        Identify the midpoints of every line segment in the triangulation.
        If an array of segments of shape (no_of_segments,2) is given,
        then the midpoints of only those segments is returned. Note,
        segments in the array must not be duplicates or the re-triangulation
        will fail. Take care not to miss that (n1,n2) is equivalent to (n2,n1).
        """

        if type(segments) == type(None):
            segments = self.identify_segments()

        mids = (self.points[segments[:,0]] + self.points[segments[:,1]]) * 0.5
        mids /= np.sqrt(mids[:,0]**2 + mids[:,1]**2 + mids[:,2]**2 ).reshape(-1,1)

        lons, lats = xyz2lonlat(mids[:,0], mids[:,1], mids[:,2])

        return lons, lats

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

        midlls = xyz2lonlat(mids[:,0], mids[:,1], mids[:,2])

        return midlls



    def lons_map_to_wrapped(self, lon):

        lons = np.array(lon)
        lons = np.mod(lon+np.pi, 2*np.pi) - np.pi

        return lons

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


    def edge_lengths(self):
        """
        Compute the edge-lengths of each triangle in the triangulation.
        """

        simplex = self.simplices.T

        # simplex is vectors a, b, c defining the corners
        a = self.points[simplex[0]]
        b = self.points[simplex[1]]
        c = self.points[simplex[2]]

        ## dot products to obtain angles
        ab = np.arccos((a * b).sum(axis=1))
        bc = np.arccos((b * c).sum(axis=1))
        ac = np.arccos((a * c).sum(axis=1))

        ## As this is a unit sphere, angle = length so ...

        return ab, bc, ac

    def angular_separation(self, p1, p2):
        """
        Compute the angles between lon / lat points p1 and p2 given in radians.
        On the unit sphere, this also corresponds to the great circle distance.
        p1 and p2 can be numpy arrays of the same length.
        """

        xp1 = lonlat2xyz(p1)
        xp2 = lonlat2xyz(p2)


        ## dot products to obtain angles

        angles = np.arccos((xp1 * xp2).sum(axis=1))

        ## As this is a unit sphere, angle = length

        return angles



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
        return points defining a refined triangulation obtained by bisection of all edges
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


    def midpoint_refine_triangulation_by_vertices(self, vertices):
        """
        return points defining a refined triangulation obtained by bisection of all edges
        in the triangulation connected to any of the vertices in the list provided
        """

        mlons, mlats = self.segment_midpoints_by_vertices(vertices=vertices)

        lonv1 = np.concatenate((self.lons, mlons), axis=0)
        latv1 = np.concatenate((self.lats, mlats), axis=0)

        return lonv1, latv1




    def edge_refine_triangulation_by_triangles(self, triangles):
        """
        return points defining a refined triangulation obtained by bisection of all edges
        in the triangulation that are associated with the triangles in the list
        of indices provided. Note that triangles are here represented as a single
        index. The vertices of triangle i are given by self.simplices[i].
        """

        ## Note there should be no duplicates in the list of triangles
        ## but because we remove duplicates from the list of all segments,
        ## there is no pressing need to check this.

        # identify the segments

        segments = []

        for index in np.array(triangles).reshape(-1):
            tri = self.simplices[index]
            segments.append( min( tuple((tri[0], tri[1])), tuple((tri[0], tri[1]))) )
            segments.append( min( tuple((tri[1], tri[2])), tuple((tri[2], tri[1]))) )
            segments.append( min( tuple((tri[0], tri[2])), tuple((tri[2], tri[0]))) )

        segs = np.array(list(set(segments)))

        mlons, mlats = self.segment_midpoints(segs)

        lonv1 = np.concatenate((self.lons, mlons), axis=0)
        latv1 = np.concatenate((self.lats, mlats), axis=0)

        return lonv1, latv1


    def edge_refine_triangulation_by_vertices(self, vertices):
        """
        return points defining a refined triangulation obtained by bisection of all edges
        in the triangulation connected to any of the vertices in the list provided
        """

        triangles = self.identify_vertex_triangles(vertices)

        return self.edge_refine_triangulation_by_triangles(triangles)



    def centroid_refine_triangulation_by_triangles(self, triangles):
        """
        return points defining a refined triangulation obtained by bisection of all edges
        in the triangulation that are associated with the triangles in the list provided.
        Note that triangles are here represented as a single
        index. The vertices of triangle i are given by self.simplices[i].
        """

        # Remove duplicates from the list of triangles

        triangles = np.unique(np.array(triangles))

        mlons, mlats = self.face_midpoints(simplices=self.simplices[triangles])

        lonv1 = np.concatenate((self.lons, mlons), axis=0)
        latv1 = np.concatenate((self.lats, mlats), axis=0)

        return lonv1, latv1


    def centroid_refine_triangulation_by_vertices(self, vertices):
        """
        return points defining a refined triangulation obtained by bisection of all edges
        in the triangulation connected to any of the vertices in the list provided
        """

        triangles = self.identify_vertex_triangles(vertices)

        return self.centroid_refine_triangulation_by_triangles(triangles)



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


    def _build_cKDtree(self):

        try:
            import scipy.spatial
            self._cKDtree =  scipy.spatial.cKDTree(self.points)

        except:
            self._cKDtree = None


    def nearest_vertices(self, lon, lat, k=1, max_distance=2.0 ):

        if self.tree == None:
            return 0, 0

        lons = np.array(lon).reshape(-1,1)
        lats = np.array(lat).reshape(-1,1)

        xyz = np.empty((lons.shape[0],3))
        x,y,z = lonlat2xyz(lons, lats)

        xyz[:,0] = x[:].reshape(-1)
        xyz[:,1] = y[:].reshape(-1)
        xyz[:,2] = z[:].reshape(-1)

        dxyz, vertices = self._cKDtree.query(xyz, k=k, distance_upper_bound=max_distance)


        if k == 1:   # force this to be a 2D array
            vertices = np.reshape(vertices, (-1, 1))

        ## Now find the angular separation / great circle distance: dlatlon


        vertxyz = self.points[vertices].transpose(0,2,1)
        extxyz  = np.repeat(xyz, k, axis=1).reshape(vertxyz.shape)

        angles = np.arccos((extxyz * vertxyz).sum(axis=1))

        return angles, vertices




## Helper functions for the module

def lonlat2xyz(lon, lat):
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

def xyz2lonlat(x,y,z):
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


def dxyz2dlonlat(x,y,z, dfx, dfy, dfz):
    """
    Take stripack df/dx, df/dy, df/dz format and convert to
    surface gradients df/dlon, df/dlat


    Notes


    """

    xs = np.array(x)
    ys = np.array(y)
    zs = np.array(z)

    lons = np.arctan2(y, x)
    lats = np.arcsin(z)

    dfxs = np.array(dfx)
    dfys = np.array(dfy)
    dfzs = np.array(dfz)

    dlon = -dfxs * np.cos(lats) * np.sin(lons) + dfys * np.cos(lats) * np.cos(lons) # no z dependence
    dlat = -dfxs * np.sin(lats) * np.cos(lons) - dfys * np.sin(lats) * np.sin(lons) + dfzs * np.cos(lats)

    corr = np.sqrt((1.0-z**2))
    valid = ~np.isclose(corr,0.0)

    dlon[valid] = dlon[valid] / corr[valid]

    return dlon , dlat

def great_circle_points(start, finish, segments):
    """
    Return a set of points along a great circle between the given points.
    start(lon, lat), finish(lon,lat) in radians.
    Fails if start and end are on a diameter
    """

    xyz = np.array(lonlat2xyz(np.array((start[0], finish[0])), np.array((start[1], finish[1]))))

    w = (np.linspace(1.0,0.0,segments))
    vectors_xyz = np.empty((segments,3))

    vectors_xyz[:,0] = xyz.T[0,0] * w[:] + xyz.T[1,0] * (1.0-w[:])
    vectors_xyz[:,1] = xyz.T[0,1] * w[:] + xyz.T[1,1] * (1.0-w[:])
    vectors_xyz[:,2] = xyz.T[0,2] * w[:] + xyz.T[1,2] * (1.0-w[:])

    vectors_xyz[:,:] /= np.sqrt((vectors_xyz[:,:]**2).sum(axis=1)).reshape(-1,1)

    lons, lats = xyz2lonlat(vectors_xyz[:,0], vectors_xyz[:,1],vectors_xyz[:,2])

    return lons, lats
