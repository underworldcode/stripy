---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

## Example 4 - `stripy` gradients on the sphere

SSRFPACK is a Fortran 77 software package that constructs a smooth interpolatory or approximating surface to data values associated with arbitrarily distributed points on the surface of a sphere. It employs automatically selected tension factors to preserve shape properties of the data and avoid overshoot and undershoot associated with steep gradients.



#### Notebook contents

   - [Analytic function and derivatives](#Analytic-function)
   - [Evaluating accuracy](#Derivatives-of-solution-compared-to-analytic-values)


The next example is [Ex5-Smoothing](./Ex5-Smoothing.md)

+++

### Define a computational mesh

Use the (usual) icosahedron with face points included.

```{code-cell} ipython3
import stripy as stripy

mesh = stripy.spherical_meshes.icosahedral_mesh(refinement_levels=4, include_face_points=True)

print(mesh.npoints)
```

### Analytic function 

Define a relatively smooth function that we can interpolate from the coarse mesh to the fine mesh and analyse

```{code-cell} ipython3
import numpy as np

def analytic(lons, lats, k1, k2):
     return np.cos(k1*lons) * np.sin(k2*lats)

def analytic_ddlon(lons, lats, k1, k2):
     return -k1 * np.sin(k1*lons) * np.sin(k2*lats) / np.cos(lats)

def analytic_ddlat(lons, lats, k1, k2):
     return k2 * np.cos(k1*lons) * np.cos(k2*lats) 

analytic_sol = analytic(mesh.lons, mesh.lats, 5.0, 2.0)
analytic_sol_ddlon = analytic_ddlon(mesh.lons, mesh.lats, 5.0, 2.0)
analytic_sol_ddlat = analytic_ddlat(mesh.lons, mesh.lats, 5.0, 2.0)
```

```{code-cell} ipython3
%matplotlib inline

import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt


fig = plt.figure(figsize=(10, 10), facecolor="none")
ax  = plt.subplot(111, projection=ccrs.Orthographic(central_longitude=0.0, central_latitude=0.0, globe=None))
ax.coastlines(color="lightgrey")
ax.set_global()

lons0 = np.degrees(mesh.lons)
lats0 = np.degrees(mesh.lats)

ax.scatter(lons0, lats0, 
            marker="o", s=10.0, transform=ccrs.PlateCarree(), c=analytic_sol, cmap=plt.cm.RdBu)

pass
```

### Derivatives of solution compared to analytic values

The `gradient_lonlat` method of the sTriangulation takes a data array reprenting values on the mesh vertices and returns the lon and lat derivatives. There is an equivalent `gradient_xyz` method which returns the raw derivatives in Cartesian form. Although this is generally less useful, if you are computing the slope (for example) that can be computed in either coordinate system and may cross the pole, consider using the Cartesian form.

```{code-cell} ipython3
stripy_ddlon, stripy_ddlat = mesh.gradient_lonlat(analytic_sol)
```

```{code-cell} ipython3
:tags: []

import k3d
plot = k3d.plot(camera_auto_fit=False, grid_visible=False, 
                menu_visibility=True, axes_helper=False )

indices = mesh.simplices.astype(np.uint32)
points = np.column_stack(mesh.points.T).astype(np.float32)

mesh_viewer = k3d.mesh(points, indices, wireframe=False, attribute=analytic_sol,
                   color_map=k3d.colormaps.basic_color_maps.CoolWarm, 
                   name="original",
                   flat_shading=False, opacity=1.0  )

plot   += mesh_viewer
plot   += k3d.points(points, point_size=0.01,color=0x779977)


plot.display()

## ## ## 

from ipywidgets import interact, interactive
import ipywidgets as widgets

choices = { "analytic": analytic_sol,
             "stripy ddlon": stripy_ddlon, 
             "stripy ddlat": stripy_ddlat, 
             "error ddlon":  stripy_ddlon-analytic_sol_ddlon, 
             "error ddlat":  stripy_ddlat-analytic_sol_ddlat }

@interact(choice=choices.keys())
def chooser(choice):
    mesh_viewer.attribute = choices[choice].astype(np.float32)
    range = np.sqrt((choices[choice]**2).mean()) * 0.5
    mesh_viewer.color_range = [-range, range]
    return 






```

The next example is [Ex5-Smoothing](./Ex5-Smoothing.md)
