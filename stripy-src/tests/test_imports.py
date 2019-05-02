from __future__ import print_function

## This would be a good test to run after a pip installation from source

## This obviously needs to be shortened but first I'm
## just adding the test in a stream-of-consciousness fashion
##
## Many of these dependencies are optional from the notebooks and
## are not part of the requirements for stripy ... they are all in
## the docker image though

print("STRIPY test_imports.py ")

## ==========================

print("\tDependency: numpy")

try:
    import numpy
    print("\t\t You have numpy version {}".format(numpy.__version__))

except:
    print("\t\t numpy import failed")


print("\tDependency: scipy")

try:
    import scipy
    print("\t\t You have scipy version {}".format(scipy.__version__))

except:
    print("\t\t scipy import failed")

## ==========================


print("\tStripy modules: stripy")

try:
    import stripy
    print("\t\t Success")

except:
    print("\t\t stripy import failed")


print("\tStripy modules: stripy.documentation")

try:
    from stripy import documentation
    print("\t\t Success")

except:
    print("\t\t stripy.documentation import failed")


print("\tStripy modules: stripy.spherical_meshes")

try:
    from stripy import spherical_meshes
    print("\t\t Success")

except:
    print("\t\t stripy.spherical_meshes import failed")


print("\tStripy modules: stripy.cartesian_meshes")

try:
    from stripy import cartesian_meshes
    print("\t\t Success")

except:
    print("\t\t stripy.cartesian_meshes import failed")



print("\tStripy modules: stripy.sTriangulation")

try:
    from stripy import sTriangulation
    print("\t\t Success")

except:
    print("\t\t stripy.sTriangulation import failed")


print("\tStripy modules: stripy.Triangulation")

try:
    from stripy import Triangulation
    print("\t\t Success")

except:
    print("\t\t stripy.Triangulation import failed")



## ============================================

print("\tJupyter notebook system")

from subprocess import check_output

try:
    result = str(check_output(['which', 'jupyter']))[2:-3]
    print("\t\t jupyter is located at {}".format(result))

    result = str(check_output(['jupyter', 'notebook', "--version"]))[2:-3]
    print("\t\t jupyter notebook version {}".format(result))

except:
    print("\t\t Unable to find jupyter notebook installation")


## ============================================

print("\tNotebook Helper: matplotlib")

try:
    import matplotlib
    print("\t\t Success")

except:
    print("\t\t matplotlib import failed")
    print("\t\t You may not be able to plot graphs in the notebooks")



print("\tNotebook Helper: cartopy")

try:
    import cartopy
    print("\t\t Success")

except:
    print("\t\t cartopy import failed")
    print("\t\t You may not be able to plot maps in the notebooks")

print("\tNotebook Helper: gdal")

try:
    import gdal
    print("\t\t Success")

except:
    print("\t\t gdal import failed")
    print("\t\t You may not be able to manipulate images for mapping in the notebooks")


print("\tNotebook Helper: lavavu")

try:
    import lavavu
    print("\t\t Success")

except:
    print("\t\t lavavu import failed")
    print("\t\t You may not be able to visualise some of the 3D mesh examples ")


print("\tNotebook Helper: pyproj")

try:
    import lavavu
    print("\t\t Success")

except:
    print("\t\t pyproj import failed")
    print("\t\t Some of the mapping examples may not work as expected ")


## ==================================


print("\tStripy use case: litho1pt0")

try:
    import litho1pt0
    print("\t\t Success")

except:
    print("\t\t litho1pt0 import failed")
    print("\t\t You may not be able to run the use case notebooks ")
