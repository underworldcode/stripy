import pytest
import stripy
import numpy as np

def test_nearest_nd_interpolation():

    p0 = 0.0
    p2 = np.pi/2
    p4 = np.pi/4
    
    coords = np.array([[p0 , -p2], \
                       [-p2,  p0], \
                       [p0 ,  p2], \
                       [p2 ,  p0]])

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats)

    Z = np.linspace(-p2, p2, mesh.npoints)

    ilons = lons[0] + 0.001
    ilats = lats[0] + 0.001

    Zi, ierr = mesh.interpolate_nearest(ilons, ilats, Z)
    
    # check if we return the nearest numbers
    if Zi == Z[0]:
        print("PASS! (Interpolation - nearest neighbour)")
    else:
        assert False, "FAIL! (Interpolation - nearest neighbour)"


def test_linear_interpolation():

    p0 = 0.0
    p2 = np.pi/2
    p4 = np.pi/4

    coords = np.array([[p0 , -p2], \
                       [-p2,  p0], \
                       [p0 ,  p2], \
                       [p2 ,  p0], \
                       [p0 ,  p0]])

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats)

    Z = mesh.lons

    npts = 5
    ilons = np.linspace(-p2, p2, npts)
    ilats = np.zeros(npts)

    Zi, ierr = mesh.interpolate_linear(ilons, ilats, Z)

    # this should be true
    # but machine precision may differ so we don't test it
    # print((Zi == ilons).all())

    bounded = Zi[0] == ilons[0] and Zi[-1] == ilons[-1]
    ascending = ( np.diff(Zi) > 0 ).all()

    # check if linear gradient across the equator
    if bounded and ascending:
        print("PASS! (Interpolation - linear")
    else:
        assert False, "FAIL! (Interpolation - linear)"


def test_cubic_interpolation():

    p0 = 0.0
    p2 = np.pi/2
    p4 = np.pi/4

    # we need more points for cubic interpolation
    coords = np.array([[p0 , -p2], \
                       [-p2,  p0], \
                       [p0 ,  p2], \
                       [p2 ,  p0], \
                       [p0 , -p4], \
                       [-p4,  p0], \
                       [p0 ,  p4], \
                       [p4 ,  p0]])

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats)

    Z = mesh.lons**2

    npts = 7
    ilons = np.linspace(-p2, p2, npts)
    ilats = np.zeros(npts)

    Zi_linear, ierr = mesh.interpolate_linear(ilons, ilats, Z)
    Zi_cubic,  ierr = mesh.interpolate_cubic(ilons, ilats, Z)

    diff_linear = np.abs(Zi_linear - ilons**2).sum()
    diff_cubic  = np.abs(Zi_cubic  - ilons**2).sum()

    # check if cubic interpolation is more accurate than linear
    if diff_cubic < diff_linear:
        print("PASS! (Interpolation - cubic")
    else:
        assert False, "FAIL! (Interpolation - cubic)"


def test_derivative():

    p0 = 0.0
    p2 = np.pi/2
    p4 = np.pi/4

    coords = np.array([[p0 , -p2], \
                       [-p2,  p0], \
                       [p0 ,  p2], \
                       [p2 ,  p0], \
                       [p0 , -p4], \
                       [-p4,  p0], \
                       [p0 ,  p4], \
                       [p4 ,  p0], \
                       [-p4, -p4], \
                       [-p4,  p4], \
                       [p4 ,  p4], \
                       [p4 , -p4], \
                       [p0 ,  p0]])
    
    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats)

    # create a soup bowl
    Z = mesh.lons**2 + mesh.lats**2

    # derivatives will have a constant gradient
    dZdlon, dZdlat = mesh.gradient_lonlat(Z, nit=10, tol=1e-12)

    # interpolate onto a straight line
    ipts = np.linspace(-p2, p2, 5)
    dZdlon_interp, ierr = mesh.interpolate_linear(ipts, ipts*0, dZdlon)
    dZdlat_interp, ierr = mesh.interpolate_linear(ipts*0, ipts, dZdlat)

    ascending_lon = ( np.diff(dZdlon_interp) > 0 ).all()
    ascending_lat = ( np.diff(dZdlat_interp) > 0 ).all()

    if ascending_lon and ascending_lat:
        print("PASS! (Derivatives)")
    else:
        assert False, "FAIL! (Derivatives)"


def test_smoothing():

    p0 = 0.0
    p2 = np.pi/2
    p4 = np.pi/4

    # we need more points for cubic interpolation
    coords = np.array([[p0 , -p2], \
                       [-p2,  p0], \
                       [p0 ,  p2], \
                       [p2 ,  p0], \
                       [p0 , -p4], \
                       [-p4,  p0], \
                       [p0 ,  p4], \
                       [p4 ,  p0], \
                       [p0 ,  p0]])

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats)

    Z = np.ones_like(lons)
    Z[-1] = 0.0

    weights = np.ones_like(lons)
    f_smooth, ierr = mesh.smoothing(Z, weights, 0.05, 1e-2, 1e-5)

    # check if f_smooth is smoother than f
    if (f_smooth.max() - f_smooth.min()) < 1.0:
        print("PASS! (Smoothing)")
    else:
        assert False, "FAIL! (Smoothing)"


if __name__ == "__main__":
    test_derivative()
    test_nearest_nd_interpolation()
    test_linear_interpolation()
    test_cubic_interpolation()
    test_smoothing()
