---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.12
    jupytext_version: 1.6.0
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

## Properties of the cratons

Requirement is to obtain some geophysical data about the cratonic lithosphere and compare to the same property for some other geological setting.

We can use the crust 1.0 geological regionalisation to query certain parts of the litho1.0 dataset and determine how lithospheric properties relate to the geological setting.

To do so, we need to be able to:

   - Construct a uniform discretisation of the globe for sampling
   - Obtain crustal type at those points
   - Select the points that meet certain criteria in the Crust 1.0 model
   - Extract properties from Litho1.0 for those points
   - Analyse
   - Plot
   
Stripy can produce the relevant sample points (we could, perhaps, use the original sample points for the litho1.0 model, or well- distributed points from a face-included icosahedron). Litho1.0 and Crust1.0 have the interpolation tools we need.

```{code-cell} ipython3
import litho1pt0 as litho1pt0
import stripy
import numpy as np
```

```{code-cell} ipython3
:hide_input: false

even_mesh = stripy.spherical_meshes.icosahedral_mesh(include_face_points=True, refinement_levels=5)
```

```{code-cell} ipython3
number_of_mesh_points = even_mesh.npoints
print (number_of_mesh_points)
print (even_mesh.lats.shape)
```

```{code-cell} ipython3
latitudes_in_radians = even_mesh.lats
latitudes_in_degrees = np.degrees(latitudes_in_radians)

longitudes_in_radians = even_mesh.lons 
longitudes_in_degrees = np.degrees(longitudes_in_radians)%360.0 - 180.0
```

```{code-cell} ipython3
## Plot these points to check the mesh is reasonable

print (latitudes_in_degrees.min(), latitudes_in_degrees.max())
print (longitudes_in_degrees.min(), longitudes_in_degrees.max())
```

```{code-cell} ipython3
## Other version of this

# Make an empty array with the same size as the number of mesh points
crustype2 = np.empty(number_of_mesh_points)
crustype2 = np.empty_like(even_mesh.lats)

# Loop and fill the array with crust regionalisation values
for i in range(0,even_mesh.npoints):
    if i%1000 == 0:
        print("Iteration: {}".format(i))
    crustype2[i] = litho1pt0.crust_type_at(lat=latitudes_in_degrees[i], lon=longitudes_in_degrees[i])
    
    
```

```{code-cell} ipython3
## Now need depth values from litho1pt0 geophysical model 

litho1pt0.l1_layer_decode
```

```{code-cell} ipython3
l1 = litho1pt0.layer_depth(latitudes_in_degrees, longitudes_in_degrees, "LID-BOTTOM") *0.001
l2 = litho1pt0.layer_depth(latitudes_in_degrees, longitudes_in_degrees, "LID-TOP") *0.001

lthickness = (l1 - l2)

print("Litho thickness - {}, {}".format(lthickness.min(), lthickness.max()))
```

```{code-cell} ipython3
c1 = litho1pt0.layer_depth(latitudes_in_degrees, longitudes_in_degrees, "CRUST3-BOTTOM") *0.001
c2 = litho1pt0.layer_depth(latitudes_in_degrees, longitudes_in_degrees, "CRUST1-TOP") *0.001

s0 = litho1pt0.layer_depth(latitudes_in_degrees, longitudes_in_degrees, "SEDS1-TOP") *0.001

cthickness = (c1 - c2)
sthickness = (c2 - s0)

print("Crust thickness - {}, {}".format(cthickness.min(), cthickness.max()))
```

```{code-cell} ipython3
## Is that sensible ?

%matplotlib inline

import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

global_extent     = [-180.0, 180.0, -89, 89]

projection0 = ccrs.PlateCarree()
projection1 = ccrs.Orthographic(central_longitude=140.0, central_latitude=0.0, globe=None)
projection2 = ccrs.Mollweide()
projection3 = ccrs.Robinson()
base_projection = ccrs.PlateCarree()

fig = plt.figure(figsize=(12, 12), facecolor="none")
ax  = plt.subplot(111, projection=projection0)
ax.set_global()

colormap = plt.cm.get_cmap('RdYlBu_r', 10)

m = ax.scatter(longitudes_in_degrees, latitudes_in_degrees, c=sthickness, cmap=colormap,
               s=2.0, linewidth=0.0, transform=ccrs.PlateCarree())

plt.colorbar(mappable=m, orientation="horizontal", shrink=0.5)

#ax.add_feature(cartopy.feature.OCEAN, alpha=0.5, zorder=99, facecolor="#BBBBBB")
ax.coastlines(resolution="50m", zorder=100, linewidth=0.5)
```

```{code-cell} ipython3
crustype = litho1pt0.crust_type_at(lat=np.degrees(even_mesh.lats), lon=np.degrees(even_mesh.lons))
```

```{code-cell} ipython3
# litho1pt0.c1_region_descriptor in 2,3,4 are Archean terranes - this is an easy way to get those values

craton_cthickness = cthickness[np.logical_and(crustype<=4, crustype>=2)]
craton_lthickness = lthickness[np.logical_and(crustype<=4, crustype>=2)]

# litho1pt0.c1_region_descriptor in 5,6,7,8 are Proterozoic terranes - this is an easy way to get those values

proton_cthickness = cthickness[np.logical_and(crustype<=8, crustype>=5)]
proton_lthickness = lthickness[np.logical_and(crustype<=8, crustype>=5)]
```

```{code-cell} ipython3
print(craton_lthickness.mean())
print(proton_lthickness.mean())
```

```{code-cell} ipython3
for i, crustName in enumerate(litho1pt0.c1_region_descriptor):
    if "Archean" in crustName:
        print(i,"-",crustName)
```

```{code-cell} ipython3

```
