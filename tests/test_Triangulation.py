import numpy as np
import stripy

from scipy import spatial
from time import clock

try: range = xrange
except: pass
    

def test_derivative(mesh):
    from scipy import interpolate

    # Create a field to test derivatives
    x, y = mesh.x, mesh.y
    Z = np.exp(-x**2 - y**2)
    Zx = -2*x*Z
    Zy = -2*y*Z
    gradZ = np.hypot(Zx, Zy)

    # Stripy
    t = clock()
    Zx1, Zy1 = mesh.gradient(Z, nit=10, tol=1e-10)
    t1 = clock() - t
    gradZ1 = np.hypot(Zx1, Zy1)

    # Spline
    spl = interpolate.SmoothBivariateSpline(x, y, Z)
    t = clock()
    Zx2 = spl.ev(x, y, dx=1)
    Zy2 = spl.ev(x, y, dy=1)
    t2 = clock() - t
    gradZ2 = np.hypot(Zx2, Zy2)

    # Clough Tocher
    # This one is most similar to what is used in stripy
    t = clock()
    cti = interpolate.CloughTocher2DInterpolator(np.column_stack([x,y]),\
                                                 Z, tol=1e-10, maxiter=20)
    t3 = clock() - t
    Zx3 = cti.grad[:,:,0].ravel()
    Zy3 = cti.grad[:,:,1].ravel()
    gradZ3 = np.hypot(Zx3, Zy3)

    res1 = ((gradZ1 - gradZ)**2).max()
    res2 = ((gradZ2 - gradZ)**2).max()
    res3 = ((gradZ3 - gradZ)**2).max()
    print("squared error in first derivative\n  \
           - stripy = {} took {}s\n  \
           - spline = {} took {}s\n  \
           - cloughtocher = {} took {}s".format(res1, t1, res2, t2, res3, t3))


def test_interpolation(mesh):
    from scipy import interpolate

    x, y = mesh.x, mesh.y
    Z = np.exp(-x**2 - y**2)
    # We ensure interpolation points are within convex hull
    xi = np.random.uniform(0.1, 0.9, 10)
    yi = np.random.uniform(0.1, 0.9, 10)

    # Stripy
    zn1 = mesh.interpolate_nearest(xi, yi, Z)
    zl1 = mesh.interpolate_linear(xi, yi, Z)
    zc1 = mesh.interpolate_cubic(xi, yi, Z)

    # cKDTree
    tree = interpolate.NearestNDInterpolator((x,y), Z)
    zn2 = tree(xi, yi)

    # Qhull
    tri = interpolate.LinearNDInterpolator((x,y), Z, 0.0)
    zl2 = tri((xi, yi))

    # Clough Tocher
    cti = interpolate.CloughTocher2DInterpolator(np.column_stack([x,y]),\
                                                 Z, tol=1e-10, maxiter=20)
    zc2 = cti((xi, yi))
    zc2[np.isnan(zc2)] = 0.0

    # Spline
    spl = interpolate.SmoothBivariateSpline(x, y, Z)
    zc3 = spl.ev(xi,yi)

    # Radial basis function
    rbf = interpolate.Rbf(x, y, Z)
    zc4 = rbf(xi, yi)

    print("squared residual in interpolation\n  \
           - nearest neighbour = {}\n  \
           - linear = {}\n  \
           - cubic (clough-tocher) = {}\n  \
           - cubic (spline) = {}\n  \
           - cubic (rbf) = {}".format(((zn1 - zn2)**2).max(), \
                                      ((zl1 - zl2)**2).max(), \
                                      ((zc1 - zc2)**2).max(), \
                                      ((zc1 - zc3)**2).max(), \
                                      ((zc1 - zc4)**2).max(),) )


def test_smoothing(mesh):
    pass

if __name__ == "__main__":
    
    # Create some (semi) random points
    np.random.seed(0)
    x = np.random.random(500)
    y = np.random.random(500)

    # triangulation
    mesh = stripy.Triangulation(x, y)

    test_derivative(mesh)
    test_interpolation(mesh)