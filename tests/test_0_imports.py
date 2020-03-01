import pytest

## ==========================

def test_numpy_import():
    import numpy
    return

def test_scipy_import():
    import scipy
    print("\t\t You have scipy version {}".format(scipy.__version__))


def test_stripy_modules():
    import stripy
    from stripy import documentation
    from stripy import spherical_meshes
    from stripy import cartesian_meshes
    from stripy import sTriangulation
    from stripy import Triangulation


def test_jupyter_available():
    from subprocess import check_output
    try:
        result = str(check_output(['which', 'jupyter']))[2:-3]
    except:
        print("Jupyter not installed")
        print("Jupyter is needed to run the example documentation")

# def test_documentation_dependencies():
#     import matplotlib
#     import cartopy
#     import imageio
#     import lavavu
#     import pyproj
#
# def test_litho1pt0():
#     import litho1pt0
