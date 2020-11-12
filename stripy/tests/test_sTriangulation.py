import pytest
import stripy
import numpy as np

@pytest.mark.parametrize("permute", [False, True])
def test_nearest_nd_interpolation(permute):

    p0 = 0.0
    p2 = np.pi/2
    p4 = np.pi/4

    coords = np.array([[p0 , -p2+np.random.random()*1.0e-16], \
                       [-p2,  p0+np.random.random()*1.0e-16], \
                       [p0 ,  p2+np.random.random()*1.0e-16], \
                       [p2 ,  p0+np.random.random()*1.0e-16]])

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats, permute=permute)

    Z = np.linspace(-p2, p2, mesh.npoints)

    ilons = lons[0] + 0.001
    ilats = lats[0] + 0.001

    Zi, ierr = mesh.interpolate_nearest(ilons, ilats, Z)

    # check if we return the nearest numbers
    if Zi == Z[0]:
        print("PASS! (Interpolation - nearest neighbour)")
    else:
        assert False, "FAIL! (Interpolation - nearest neighbour)"


@pytest.mark.parametrize("permute", [False, True])
def test_linear_interpolation(permute):

    p0 = 0.0
    p2 = np.pi/2
    p4 = np.pi/4

    coords = np.array([[p0 , -p2+np.random.random()*1.0e-16], \
                       [-p2,  p0+np.random.random()*1.0e-16], \
                       [p0 ,  p2+np.random.random()*1.0e-16], \
                       [p2 ,  p0+np.random.random()*1.0e-16], \
                       [p0 ,  p0+np.random.random()*1.0e-16]])

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats, permute=permute)

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


@pytest.mark.parametrize("permute", [False, True])
def test_cubic_interpolation(permute):

    p0 = 0.0
    p2 = np.pi/2
    p4 = np.pi/4

    # we need more points for cubic interpolation
    coords = np.array([[p0 , -p2+np.random.random()*1.0e-16], \
                       [-p2,  p0+np.random.random()*1.0e-16], \
                       [p0 ,  p2+np.random.random()*1.0e-16], \
                       [p2 ,  p0+np.random.random()*1.0e-16], \
                       [p0 , -p4+np.random.random()*1.0e-16], \
                       [-p4,  p0+np.random.random()*1.0e-16], \
                       [p0 ,  p4+np.random.random()*1.0e-16], \
                       [p4 ,  p0+np.random.random()*1.0e-16]])

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats, permute=permute)

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


@pytest.mark.parametrize("permute", [False, True])
def test_cubic_interpolation_tension(permute):

    p0 = 0.0
    p2 = np.pi/2
    p4 = np.pi/4

    # we need more points for cubic interpolation
    coords = np.array([[p0 , -p2+np.random.random()*1.0e-16], \
                       [-p2,  p0+np.random.random()*1.0e-16], \
                       [p0 ,  p2+np.random.random()*1.0e-16], \
                       [p2 ,  p0+np.random.random()*1.0e-16], \
                       [p0 , -p4+np.random.random()*1.0e-16], \
                       [-p4,  p0+np.random.random()*1.0e-16], \
                       [p0 ,  p4+np.random.random()*1.0e-16], \
                       [p4 ,  p0+np.random.random()*1.0e-16]])

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats, permute=permute)

    Z = mesh.lons**2

    npts = 7
    ilons = np.linspace(-p2, p2, npts)
    ilats = np.zeros(npts)

    Zi_linear, ierr = mesh.interpolate_linear(ilons, ilats, Z)
    Zi_cubic,  ierr = mesh.interpolate_cubic(ilons, ilats, Z)

    sigma = mesh.get_spline_tension_factors(Z)
    Zi_cubicT,  ierr = mesh.interpolate_cubic(ilons, ilats, Z, sigma=sigma)

    sigma.fill(45.)
    Zi_cubicTmax,  ierr = mesh.interpolate_cubic(ilons, ilats, Z, sigma=sigma)


    if np.abs(Zi_cubicT - Zi_cubic).any():
        print("PASS! (Interpolation - cubic tensioned splines")
    # check if cubic interpolation with max tension is like linear interpolation
    elif np.abs(Zi_linear-Zi_cubicTmax).sum() < np.abs(Zi_cubic-Zi_cubicTmax).sum():
        print("PASS! (Interpolation - cubic tensioned splines")
    else:
        assert False, "FAIL! (Interpolation - cubic tensioned splines)"


@pytest.mark.parametrize("permute", [False])
def test_cubic_interpolation_grid(permute):

    p0 = 0.0
    p2 = np.pi/2
    p4 = np.pi/4

    # we need more points for cubic interpolation
    coords = np.array([[p0 , -p2+np.random.random()*1.0e-16], \
                       [-p2,  p0+np.random.random()*1.0e-16], \
                       [p0 ,  p2+np.random.random()*1.0e-16], \
                       [p2 ,  p0+np.random.random()*1.0e-16], \
                       [p0 , -p4+np.random.random()*1.0e-16], \
                       [-p4,  p0+np.random.random()*1.0e-16], \
                       [p0 ,  p4+np.random.random()*1.0e-16], \
                       [p4 ,  p0+np.random.random()*1.0e-16]])

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats, permute=permute)

    Z = mesh.lons**2

    npts = 7
    ilons = np.linspace(-p2, p2, npts)
    ilats = np.linspace(-p2, p2, npts)
    lonq, latq = np.meshgrid(ilons, ilats)
    shape = (npts, npts)

    Zi_cubic, zierr = mesh.interpolate_cubic(lonq.ravel(), latq.ravel(), Z)
    Zi_cubic_grid   = mesh.interpolate_to_grid(ilons, ilats, Z)

    sigma = mesh.get_spline_tension_factors(Z, tol=1e-5)
    Zi_cubic_grid_S = mesh.interpolate_to_grid(ilons, ilats, Z, sigma=sigma)

    err_msg = "Interpolate to grid - cubic tensioned splines"

    np.testing.assert_allclose(Zi_cubic.reshape(shape), Zi_cubic_grid, atol=0.1, err_msg=err_msg)
    np.testing.assert_allclose(Zi_cubic_grid_S, Zi_cubic_grid, atol=0.5, err_msg=err_msg)

    # This has an unpredictable error in it that is causing trouble
    # for all of our automated tests and conda builds:

    assert (Zi_cubic_grid_S != Zi_cubic_grid).any(), err_msg

    # unstructured and grid interpolation works
    # and applying tension alters the result


@pytest.mark.parametrize("permute", [False, True])
def test_derivative(permute):

    p0 = 0.0
    p2 = np.pi/2
    p4 = np.pi/4

    coords = np.array([[p0 , -p2+np.random.random()*1.0e-16], 
                       [-p2,  p0+np.random.random()*1.0e-16], 
                       [p0 ,  p2+np.random.random()*1.0e-16], 
                       [p2 ,  p0+np.random.random()*1.0e-16], 
                       [p0 , -p4+np.random.random()*1.0e-16], 
                       [-p4,  p0+np.random.random()*1.0e-16], 
                       [p0 ,  p4+np.random.random()*1.0e-16], 
                       [p4 ,  p0+np.random.random()*1.0e-16], 
                       [-p4, -p4+np.random.random()*1.0e-16], 
                       [-p4,  p4+np.random.random()*1.0e-16], 
                       [p4 ,  p4+np.random.random()*1.0e-16], 
                       [p4 , -p4+np.random.random()*1.0e-16], 
                       [p0 ,  p0+np.random.random()*1.0e-16]])

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats, permute=permute)

    # create a soup bowl
    Z = mesh.lons**2 + mesh.lats**2

    # derivatives will have a constant gradient
    dZdlon, dZdlat = mesh.gradient_lonlat(Z, nit=10, tol=1e-12)

    # interpolate onto a straight line
    ipts = np.linspace(-p2, p2, 5)
    dZdlon_interp, ierr = mesh.interpolate_linear(ipts, ipts*0, dZdlon)
    dZdlat_interp, ierr = mesh.interpolate_linear(ipts*0, ipts, dZdlat)

    ascending_lon = ( np.diff(dZdlon_interp) > 1.0e-16 ).all()
    ascending_lat = ( np.diff(dZdlat_interp) > 1.0e-16 ).all()

    if ascending_lon and ascending_lat:
        print("PASS! (Derivatives)")
    else:
        assert False, "FAIL! (Derivatives)"


@pytest.mark.parametrize("permute", [False, True])
def test_smoothing(permute):

    p0 = 0.0
    p2 = np.pi/2
    p4 = np.pi/4

    # we need more points for cubic interpolation
    coords = np.array([[p0 , -p2+np.random.random()*1.0e-16], \
                       [-p2,  p0+np.random.random()*1.0e-16], \
                       [p0 ,  p2+np.random.random()*1.0e-16], \
                       [p2 ,  p0+np.random.random()*1.0e-16], \
                       [p0 , -p4+np.random.random()*1.0e-16], \
                       [-p4,  p0+np.random.random()*1.0e-16], \
                       [p0 ,  p4+np.random.random()*1.0e-16], \
                       [p4 ,  p0+np.random.random()*1.0e-16]])

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats, permute=permute)

    Z = np.ones_like(lons)
    Z[-1] = 0.0

    weights = np.ones_like(lons)
    f_smooth, f_derivatives, ierr = mesh.smoothing(Z, weights, 0.05, 1e-2, 1e-5)

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
    test_cubic_interpolation_tension()
    test_cubic_interpolation_grid()
    test_smoothing()
