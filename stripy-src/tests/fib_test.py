from __future__ import print_function
try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap
    has_matplotlib = True
except:
    has_matplotlib = False
import numpy as np
from stripack import trmesh
import time

def fibonacci_pts(npts):
    # return lats and lons of N=npts fibonacci grid on a sphere.
    pi = np.pi
    inc = pi * (3.0 - np.sqrt(5.0))
    off = 2. / npts
    lats = []; lons = []
    for k in range(npts):
       y = k*off - 1. + 0.5*off
       r = np.sqrt(1 - y**2)
       phi = k * inc
       x = np.cos(phi)*r
       z = np.sin(phi)*r
       theta = np.arctan2(np.sqrt(x**2+y**2),z)
       phi = np.arctan2(y,x)
       lats.append( 0.5*pi-theta )
       if phi < 0.: phi = 2.*pi+phi
       lons.append( phi )
    return np.array(lats), np.array(lons)


# fake test data.
def test_func(lon, lat):
    nexp = 8
    return np.cos(nexp*lon)*np.sin(0.5*lon)**nexp*np.cos(lat)**nexp+np.sin(lat)**nexp
npts = 100000
lats, lons = fibonacci_pts(npts)
shuffle = True # shuffling the points speeds up the triangulation
if shuffle:
    ix = np.arange(npts)
    np.random.shuffle(ix)
    lats = lats[ix]; lons = lons[ix]
icos_data = test_func(lons,lats)

t1 = time.clock()
print('triangulation of', len(lons),' points')
tri = trmesh(lons, lats)
print('triangulation took',time.clock()-t1,' secs')

nlons = 1440; nlats = nlons/2 + 1 # 1 degree output mesh
delta = 360./nlons
olons = delta*np.arange(nlons)
olats = -90.0 + delta*np.arange(nlats)
olons = np.radians(olons);  olats = np.radians(olats)
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
