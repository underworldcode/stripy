import pytest
import stripy
import numpy as np

def test_derivative():
    def analytic(lons, lats, k1, k2):
        return np.cos(k1*lons) * np.sin(k2*lats)

    def analytic_ddlon(lons, lats, k1, k2):
        return -k1 * np.sin(k1*lons) * np.sin(k2*lats) / np.cos(lats)

    def analytic_ddlat(lons, lats, k1, k2):
        return k2 * np.cos(k1*lons) * np.cos(k2*lats)


    mesh = stripy.spherical_meshes.icosahedral_mesh(include_face_points=True, refinement_levels=3)
    
    Z = analytic(mesh.lons, mesh.lats, 5.0, 2.0)
    dZdlon, dZdlat = mesh.gradient_lonlat(Z, nit=10, tol=1e-12)

    dZdlon_analytic = analytic_ddlon(mesh.lons, mesh.lats, 5.0, 2.0)
    dZdlat_analytic = analytic_ddlat(mesh.lons, mesh.lats, 5.0, 2.0)    

    if np.isclose(dZdlon, dZdlon_analytic, 1.,1.).all() and \
       np.isclose(dZdlat, dZdlat_analytic, 1.,1.).all():
        print("PASS! (Derivatives)")
    else:
        assert False, "FAIL! (Derivatives)"


def test_nearest_nd_interpolation():
    
    coords = np.array([[0.0, 0.0], \
                       [0.0, 1.0], \
                       [1.0, 0.0], \
                       [1.0, 1.0]])

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats)

    Z = mesh.lons + mesh.lats

    Zi, ierr = mesh.interpolate_nearest(0.9, 0.9, Z)
    
    if Zi == 2.0:
        print("PASS! (Interpolation - nearest neighbour)")
    else:
        assert False, "FAIL! (Interpolation - nearest neighbour)"


def test_linear_interpolation():

    coords = np.array([[0.0, 0.0], \
                       [0.0, 1.0], \
                       [1.0, 0.0], \
                       [1.0, 1.0]])

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats)

    Z = mesh.lons + mesh.lats

    Zi, ierr = mesh.interpolate_linear(0.5, 0.5, Z)

    if np.isclose(Zi, 1.05527):
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

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats)

    Z = mesh.lons + mesh.lats

    Zi, ierr = mesh.interpolate_cubic(0.5, 0.5, Z)

    if np.isclose(Zi, 0.94596):
        print("PASS! (Interpolation - cubic")
    else:
        assert False, "FAIL! (Interpolation - cubic)"


def test_smoothing():

    coords = np.array([[0.0, 0.0], \
                       [0.0, 1.0], \
                       [1.0, 0.0], \
                       [1.0, 1.0], \
                       [0.5, 0.5]])

    lons, lats = coords[:,0], coords[:,1]
    mesh = stripy.sTriangulation(lons, lats)

    Z = np.ones_like(lons)
    Z[-1] = 0.0

    weights = np.ones_like(lons)
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
