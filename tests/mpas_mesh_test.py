from __future__ import print_function
try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    has_matplotlib = True
except:
    has_matplotlib = False
from netCDF4 import Dataset
import numpy as np
from stripack import trmesh
import time

# MPAS model mesh files from
# https://mpas-dev.github.io/atmosphere/atmosphere_meshes.html
mesh_filename = 'x1.2562.grid.nc' # 480 km mesh
#mesh_filename = 'x1.40962.grid.nc' # 120 km mesh
#mesh_filename = 'x1.2621442.grid.nc' # 15 km mesh
mesh_nc = Dataset(mesh_filename)
lats = mesh_nc.variables['latCell'][:]
print('min/max lats:',lats.min(), lats.max())
lons = mesh_nc.variables['lonCell'][:]
print('min/max lons:',lons.min(), lons.max())

# fake test data.
def test_func(lon, lat):
    nexp = 8
    return np.cos(nexp*lon)*np.sin(0.5*lon)**nexp*np.cos(lat)**nexp+np.sin(lat)**nexp
icos_data = test_func(lons,lats)

t1 = time.clock()
print('triangulation of', len(lons),' points')
tri = trmesh(lons, lats)
print('triangulation took',time.clock()-t1,' secs')

nlons = 360; nlats = nlons/2 # 1 degree output mesh
delta = 360./nlons
olons = delta*np.arange(nlons)
olats = -90.0 + 0.5*delta + delta*np.arange(nlats)
olons = np.radians(olons)
olats = np.radians(olats)
olons, olats = np.meshgrid(olons, olats)

t1 = time.clock()
order = 1 # can be 0 (nearest neighbor) or 1 (linear)
latlon_data = tri.interp(olons,olats,icos_data,order=order)
print('interpolation took',time.clock()-t1,' secs')

latlon_datax = test_func(olons,olats)
print('max abs error:',(np.abs(latlon_datax-latlon_data)).max())
print('min/max field:',latlon_data.min(), latlon_datax.max())

# make plot on output mesh
if has_matplotlib:
    fig = plt.figure(figsize=(12,6))
    fig.add_subplot(1,2,1)
    m = Basemap(projection='ortho',lon_0=180,lat_0=20)
    x,y = m(np.degrees(olons), np.degrees(olats))
    m.drawcoastlines()
    m.drawmapboundary()
    m.contourf(x,y,latlon_data,15)
    plt.title('interpolated field order=%s' % order)
    m.colorbar()
    fig.add_subplot(1,2,2)
    m.drawcoastlines()
    m.drawmapboundary()
    m.contourf(x,y,latlon_data-latlon_datax,15)
    plt.title('error')
    m.colorbar()
    plt.tight_layout()
    plt.show()
