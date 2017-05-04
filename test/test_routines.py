import numpy as np
import stripy
from time import clock


def time_routine(routine, *args):
	t = clock()
	routine(*args)
	str_fmt = "- {:25} completed in  {:1.6f} seconds"
	print(str_fmt.format(str(routine.__name__), clock()-t))

"""
Test Triangulation routines
"""
print("---\nTriangulation routines\n---")

# Create some (semi) random points
np.random.seed(0)
x = np.random.random(10000)
y = np.random.random(10000)

z = np.hypot(x,y)

xi = np.random.random(100)
yi = np.random.random(100)

# triangulation
tri = stripy.Triangulation(x, y)

time_routine(tri.__init__, x, y)
time_routine(tri.areas)
time_routine(tri.convex_hull)
time_routine(tri.gradient, z)
time_routine(tri.gradient_local, z, 0)
time_routine(tri.interpolate_nearest, xi, yi, z)
time_routine(tri.interpolate_linear, xi, yi, z)
time_routine(tri.interpolate_cubic, xi, yi, z)
time_routine(tri.nearest_neighbour, 0.5, 0.5)
time_routine(tri.neighbour_simplices)
time_routine(tri.neighbour_and_arc_simplices)


"""
Test  sTriangulation routines
"""
print("---\nsTriangulation routines\n---")

# Create some (semi) random points
np.random.seed(0)
lons = 2.*np.pi*np.random.random(10000)
lats = np.arccos(2.*np.random.random(10000) - 1.) - np.pi/2
Z = np.exp(-lons**2 - lats**2)

lon = 2.*np.pi*np.random.random(100)
lat = np.arccos(2.*np.random.random(100) - 1.) - np.pi/2

stri = stripy.sTriangulation(lons, lats)

time_routine(stri.__init__, lons, lats)
time_routine(stri.areas)
time_routine(stri.gradient, Z)
time_routine(stri.interpolate_nearest, lon, lat, Z)
time_routine(stri.interpolate_linear, lon, lat, Z)
time_routine(stri.interpolate_cubic, lon, lat, Z)
time_routine(stri.transform_to_spherical, Z/stri.x, Z/stri.y, Z/stri.z)
time_routine(stri.tri_area, lon[:3], lat[:3])
