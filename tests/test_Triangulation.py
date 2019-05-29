import pytest
import stripy
import numpy as np

def test_derivative():
    def analytic(x, y, k1, k2):
        return np.cos(k1*x) * np.sin(k2*y)

    def analytic_dx(x, y, k1, k2):
        return -k1 * np.sin(k1*x) * np.sin(k2*y)

    def analytic_dy(x, y, k1, k2):
        return k2 * np.cos(k1*x) * np.cos(k2*y)


    extent = [0.0, 1.0, 0.0, 1.0]
    mesh = stripy.cartesian_meshes.square_mesh(extent, spacingX=0.05, spacingY=0.05, refinement_levels=1)
    
    Z = analytic(mesh.x, mesh.y, 5.0, 2.0)
    dZx, dZy = mesh.gradient(Z, nit=10, tol=1e-12)

    dZx_analytic = analytic_dx(mesh.x, mesh.y, 5.0, 2.0)
    dZy_analytic = analytic_dy(mesh.x, mesh.y, 5.0, 2.0)    

    if np.isclose(dZx, dZx_analytic, 1.,1.).all() and \
       np.isclose(dZy, dZy_analytic, 1.,1.).all():
        print("PASS! (Derivatives)")
    else:
        assert False, "FAIL! (Derivatives)"


def test_nearest_nd_interpolation():
    
    coords = np.array([[0.0, 0.0], \
                       [0.0, 1.0], \
                       [1.0, 0.0], \
                       [1.0, 1.0]])

    x, y = coords[:,0], coords[:,1]
    mesh = stripy.Triangulation(x, y)

    Z = mesh.x + mesh.y

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

    x, y = coords[:,0], coords[:,1]
    mesh = stripy.Triangulation(x, y)

    Z = mesh.x + mesh.y

    Zi, ierr = mesh.interpolate_linear(0.5, 0.5, Z)

    if np.isclose(Zi, 1.0, 0.001):
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

    Z = mesh.x + mesh.y

    Zi, ierr = mesh.interpolate_cubic(0.5, 0.5, Z)

    if np.isclose(Zi, 1.0, 0.001):
        print("PASS! (Interpolation - cubic")
    else:
        assert False, "FAIL! (Interpolation - cubic)"


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
