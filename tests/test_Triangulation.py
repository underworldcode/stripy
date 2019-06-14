import pytest
import stripy
import numpy as np

def test_nearest_nd_interpolation():
    
    coords = np.array([[0.0, 0.0], \
                       [0.0, 1.0], \
                       [1.0, 0.0], \
                       [1.0, 1.0]])

    x, y = coords[:,0], coords[:,1]
    mesh = stripy.Triangulation(x, y)

    Z = np.linspace(0.0, 10.0, mesh.npoints)

    ix = x[0] + 0.001
    iy = y[0] + 0.001

    Zi, ierr = mesh.interpolate_nearest(ix, iy, Z)
    
    if Zi == Z[0]:
        print("PASS! (Interpolation - nearest neighbour)")
    else:
        assert False, "FAIL! (Interpolation - nearest neighbour)"


def test_linear_interpolation():

    coords = np.array([[0.0, 0.0], \
                       [0.0, 1.0], \
                       [1.0, 0.0], \
                       [1.0, 1.0], \
                       [0.5, 0.5]])

    x, y = coords[:,0], coords[:,1]
    mesh = stripy.Triangulation(x, y)

    Z = mesh.x

    npts = 5
    ix = np.linspace(0.0, 1.0, npts)
    iy = np.zeros(npts)

    Zi, ierr = mesh.interpolate_linear(ix, iy, Z)

    # this should be true
    # but machine precision may differ so we don't test it
    # print((Zi == ix).all())

    bounded = Zi[0] == ix[0] and Zi[-1] == ix[-1]
    ascending = ( np.diff(Zi) > 0 ).all()

    if bounded and ascending:
        print("PASS! (Interpolation - linear")
    else:
        assert False, "FAIL! (Interpolation - linear)"


def test_cubic_interpolation():

    # we need more points for cubic interpolation
    coords = np.array([[0.0, 0.0], \
                       [0.0, 1.0], \
                       [1.0, 0.0], \
                       [1.0, 1.0], \
                       [0.1, 0.1], \
                       [0.1, 0.9], \
                       [0.9, 0.1], \
                       [0.9, 0.9]])

    x, y = coords[:,0], coords[:,1]
    mesh = stripy.Triangulation(x, y)

    Z = mesh.x**2

    npts = 7
    ilons = np.linspace(0.0, 1.0, npts)
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
    p1 = 1.0
    p2 = 2.0

    coords = np.array([[p0 , -p2], \
                       [-p2,  p0], \
                       [p0 ,  p2], \
                       [p2 ,  p0], \
                       [p0 , -p1], \
                       [-p1,  p0], \
                       [p0 ,  p1], \
                       [p1 ,  p0], \
                       [-p1, -p1], \
                       [-p1,  p1], \
                       [p1 ,  p1], \
                       [p1 , -p1], \
                       [p0 ,  p0]])
    
    x, y = coords[:,0], coords[:,1]
    mesh = stripy.Triangulation(x, y)

    # create a soup bowl
    Z = mesh.x**2 + mesh.y**2

    # derivatives will have a constant gradient
    dZdx, dZdy = mesh.gradient(Z, nit=10, tol=1e-12)

    # interpolate onto a straight line
    ipts = np.linspace(-p2, p2, 5)
    dZdx_interp, ierr = mesh.interpolate_linear(ipts, ipts*0, dZdx)
    dZdy_interp, ierr = mesh.interpolate_linear(ipts*0, ipts, dZdy)

    ascending_lon = ( np.diff(dZdx_interp) > 0 ).all()
    ascending_lat = ( np.diff(dZdy_interp) > 0 ).all()

    if ascending_lon and ascending_lat:
        print("PASS! (Derivatives)")
    else:
        assert False, "FAIL! (Derivatives)"


def test_smoothing():

    coords = np.array([[0.0, 0.0], \
                       [0.0, 1.0], \
                       [1.0, 0.0], \
                       [1.0, 1.0], \
                       [0.5, 0.5]])

    x, y = coords[:,0], coords[:,1]
    mesh = stripy.Triangulation(x, y)

    Z = np.ones_like(x)
    Z[-1] = 0.0

    weights = np.ones_like(x)
    f_smooth, ierr = mesh.smoothing(Z, weights, 0.05, 1e-2, 1e-5)

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
