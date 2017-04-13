import _stripack
import _ssrfpack
import numpy as np

try: range = xrange
except: pass

__version__ = "1.1"

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
     x : ndarray of floats, shape (n,)
        stored Cartesian x coordinates from input
     y : ndarray of floats, shape (n,)
        stored Cartesian y coordinates from input
     simplices : ndarray of ints, shape (nsimplex, 3)
        indices of the points forming the simplices in the triangulation
        points are ordered anticlockwise
     lst : ndarray of ints, shape (6n-12,)
        nodal indices with lptr and lend, define the triangulation as a set of N
        adjacency lists; counterclockwise-ordered sequences of neighboring nodes
        such that the first and last neighbors of a boundary node are boundary
        nodes (the first neighbor of an interior node is arbitrary).  In order to
        distinguish between interior and boundary nodes, the last neighbor of
        each boundary node is represented by the negative of its index.
        The indices are 1-based (as in Fortran), not zero based (as in python).
     lptr : ndarray of ints, shape (6n-12),)
        set of pointers in one-to-one correspondence with the elements of lst.
        lst(lptr(i)) indexes the node which follows lst(i) in cyclical
        counterclockwise order (the first neighbor follows the last neighbor).
        The indices are 1-based (as in Fortran), not zero based (as in python).
     lend : ndarray of ints, shape (n,)
        N pointers to adjacency lists.
        lend(k) points to the last neighbor of node K. 
        lst(lend(K)) < 0 if and only if K is a boundary node.
        The indices are 1-based (as in Fortran), not zero based (as in python).
    """
    def __init__(self, lons, lats):

        self._check_integrity(lons, lats)

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


    def _check_integrity(self, lons, lats):
        """
        Ensure lons and lats are:
         - 1D numpy arrays
         - equal size
         - within the appropriate range in radians
        """
        if len(lons.shape) != 1 or len(lats.shape) != 1:
            raise ValueError('lons and lats must be 1D')
        if lats.size != lons.size:
            raise ValueError('lons and lats must have same length')
        if (np.abs(lons)).max() > 2.*np.pi:
            raise ValueError("lons must be in radians (-2*pi <= lon <= 2*pi)")
        if (np.abs(lats)).max() > 0.5*np.pi:
            raise ValueError("lats must be in radians (-pi/2 <= lat <= pi/2)")


    def interpolate(self, lons, lats, data, order=1):

        self._check_integrity(lons, lats)

        if order not in [0,1,3]:
            raise ValueError("order must be 0, 1, or 3")
        if data.size != self.npoints:
            raise ValueError("data must be of size {}".format(self.npoints))

        interp, ierr = _stripack.interp_n(order, lats, lons,\
                                          self.x, self.y, self.z, data,\
                                          self.lst, self.lptr, self.lend)

        if ierr != 0:
            raise ValueError('ierr={} in interp_n\n{}'.format(ierr, _ier_codes[ierr]))

        return interp


    def interpolate_nearest(self,lons,lats,data):
        """
        Interpolate using nearest-neighbour approximation
        Returns the same as interpolate(lons,lats,data,order=0)
        """
        return self.interp(lons, lats, data, order=0)

    def interpolate_linear(self,lons,lats,data):
        """
        Interpolate using nearest-neighbour approximation
        Returns the same as interpolate(lons,lats,data,order=1)
        """
        return self.interp(lons, lats, data, order=1)

    def interpolate_cubic(self,lons,lats,data):
        """
        Interpolate using nearest-neighbour approximation
        Returns the same as interpolate(lons,lats,data,order=3)
        """
        return self.interp(lons, lats, data, order=3)


    def tri_area(self, tr):

        # Note the indexing must be 0 based for the py arrays
        lons = self.lons[tr-1]
        lats = self.lats[tr-1]


        xyz1 = _stripack.trans(lats[0], lons[0],1)
        xyz2 = _stripack.trans(lats[1], lons[1],1)
        xyz3 = _stripack.trans(lats[2], lons[2],1)

        area = _stripack.areas(xyz1, xyz2, xyz3)

        return area


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
        use_sigma_array = 0

        ierr = 1
        while ierr == 1:
            ierr = _ssrfpack.gradg(self.x, self.y, self.z, f, self.lst, self.lptr, self.lend,\
                                   use_sigma_array, sigma, nit, tol, gradient)
            if not guarantee_convergence:
                break

        if ierr < 0:
            raise ValueError('ierr={} in gradg\n{}'.format(ierr, _ier_codes[ierr]))

        return gradient[0], gradient[1], gradient[2]


    def find_point(self, lon, lat):

        lat = np.array(lat).astype(np.float64,copy=False)
        lon = np.array(lon).astype(np.float64,copy=False)

        npts = len(lon)

        p = np.array(_stripack.trans(lat, lon, npts)).T

        # if (np.abs(olons1)).max() > 2.*np.pi:
        #     msg="lons must be in radians (-2*pi <= lon <= 2*pi)"
        #     raise ValueError(msg)
        #
        # if (np.abs(olats1)).max() > 0.5*np.pi:
        #     msg="lats must be in radians (-pi/2 <= lat <= pi/2)"
        #     raise ValueError(msg)

        # This works one at a time only ... right now

        bccList   = []
        nodesList = []

        for i in range(0,npts):

            point = p[i]

            pt = _stripack.trfind(1,point,self.x,self.y,self.z,self.lst,self.lptr,self.lend, self.npts)

            bcc = np.array(pt[0:3], dtype=float)
            bcc /=  bcc.sum()

            # The minimum of the 3 cyclic permutations is always returned
            ttup1 = tuple((pt[3], pt[4], pt[5]))
            ttup2 = tuple((pt[4], pt[5], pt[3]))
            ttup3 = tuple((pt[5], pt[3], pt[4]))
            nodes= min(ttup1, ttup2, ttup3)

            bccList.append(bcc)
            nodesList.append(nodes)

        # if ierr != 0:
        #     raise ValueError('ierr = %s in intrpc0_n' % ierr)

        return bccList, nodesList
