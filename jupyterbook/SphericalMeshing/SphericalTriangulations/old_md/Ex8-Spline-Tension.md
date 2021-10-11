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

# Example 8 - Spline Tension

Apply spline tension to interpolaton, gradient, derivative, or smoothing routines to reduce the occurance of undershoot / overshoot inconsistencies in the solution.

The values to determine the degree of tension is stored in `sigma`. Using the routine `get_spline_tension_factors` will determine the smallest tension factor such that the spline preserves the local shape properties (monotonicity and convexity) of the data. If `sigma` is zero everywhere, then no tension is active.

We walk through a number of routines that we have explored in previous notebooks, but in this case demonstrating the use of tensioned splines.

## Contents

- [Smoothing with tension](#Smoothing-with-tension)
- [Interpolation with tension](#Interpolation-with-tension)
- [Gradients with tension](#Gradients-with-tension)

```{code-cell} ipython3
import stripy as stripy
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
%matplotlib inline
```

```{code-cell} ipython3
mesh = stripy.spherical_meshes.icosahedral_mesh(refinement_levels=4)
```

```{code-cell} ipython3
def analytic(lons, lats, k1, k2):
     return  np.cos(k1*lons) * np.sin(k2*lats) 

def analytic_noisy(lons, lats, k1, k2, noise, short):
     return  np.cos(k1*lons) * np.sin(k2*lats) + short * (np.cos(k1*5.0*lons) * np.sin(k2*5.0*lats)) +  noise * np.random.random(lons.shape)

data   = analytic(mesh.lons, mesh.lats, 5.0, 2.0)
data_n = analytic_noisy(mesh.lons, mesh.lats, 5.0, 2.0, 0.1, 0.0)
```

```{code-cell} ipython3
# get tension factors from the data
sigma   = mesh.get_spline_tension_factors(data, tol=1e-6)
sigma_n = mesh.get_spline_tension_factors(data_n, tol=1e-6)
```

## Smoothing with tension

Tension is applied simply by supplying `sigma`. It's effect is most noticible near the poles where there are edge artefacts in the solution.

```{code-cell} ipython3
stripy_smoothed,  dds, err = mesh.smoothing(data, np.ones_like(data_n), 10.0, 0.1, 0.01)
stripy_smoothed2, dds, err = mesh.smoothing(data, np.ones_like(data_n), 10.0, 0.1, 0.01, sigma=sigma)
```

```{code-cell} ipython3
import k3d

plot = k3d.plot(camera_auto_fit=False, grid_visible=False, 
                menu_visibility=True, axes_helper=False )

indices = mesh.simplices.astype(np.uint32)
points = np.column_stack(mesh.points.T).astype(np.float32)

mesh_viewer = k3d.mesh(points, indices, wireframe=False, attribute=data_n,
                   color_map=k3d.colormaps.paraview_color_maps.Cool_to_Warm, 
                   name="splines",
                   flat_shading=False, opacity=1.0  )


plot   += mesh_viewer
plot += background
plot.display()

## ## ## 

from ipywidgets import interact, interactive
import ipywidgets as widgets

choices = { "original": data_n,
             "smoothed": stripy_smoothed, 
             "smooth with tension": stripy_smoothed2, 
             "delta": stripy_smoothed - stripy_smoothed2  }

@interact(choice=choices.keys())
def chooser(choice):
    mesh_viewer.attribute = choices[choice].astype(np.float32)
    range = choices[choice].max() * 0.2
    mesh_viewer.color_range = [0, range]
    return 

```

## Interpolation with tension

Only applies to cubic interpolation. The effect of supplying a non-negative `sigma` is to produce a _more linear_ interpolation. For regions that lie outside the hull, cubic extrapolation produces wild oscillations in the solution which can be mitigated with spline tension.

```{code-cell} ipython3
# set up a discontinuous mesh

mask_points = mesh.lats < np.pi/3
cmesh = stripy.sTriangulation(mesh.lons[mask_points], mesh.lats[mask_points])

cdata   = analytic(cmesh.lons, cmesh.lats, 5.0, 2.0)
csigma  = cmesh.get_spline_tension_factors(cdata, tol=1e-6)
```

```{code-cell} ipython3
grid_z1, ierr = cmesh.interpolate_cubic(mesh.lons, mesh.lats, cdata) # no tension
grid_z2, ierr = cmesh.interpolate_cubic(mesh.lons, mesh.lats, cdata, sigma=csigma) # tension
```

```{code-cell} ipython3
import k3d

plot = k3d.plot(camera_auto_fit=False, grid_visible=False, 
                menu_visibility=True, axes_helper=False )

indices = mesh.simplices.astype(np.uint32)
points = np.column_stack(mesh.points.T).astype(np.float32)

mesh_viewer = k3d.mesh(points, indices, wireframe=False, attribute=data_n,
                   color_map=k3d.colormaps.paraview_color_maps.Cool_to_Warm, 
                   name="splines",
                   flat_shading=False, opacity=1.0  )


plot += mesh_viewer
#plot += background
plot.display()

## ## ## 

from ipywidgets import interact, interactive
import ipywidgets as widgets

choices = {  
             "cubic": grid_z1, 
             "cubic with tension": grid_z2, 
             "delta": grid_z1 - grid_z2  }

@interact(choice=choices.keys())
def chooser(choice):
    mesh_viewer.attribute = choices[choice].astype(np.float32)
    range = choices[choice].max() * 0.2
    mesh_viewer.color_range = [-range, range]
    return 


```

## Gradients with tension

Pass `sigma` to the following routines that involve derivatives:

- `gradient_lonlat`
- `gradient_xyz`
- `derivatives_lonlat`

Again, the largest difference is visible at the poles.

```{code-cell} ipython3
dlon1, dlat1 = mesh.gradient_lonlat(data, nit=5, tol=1e-6) # no tension
dlon2, dlat2 = mesh.gradient_lonlat(data, nit=5, tol=1e-6, sigma=sigma) # tension
```

```{code-cell} ipython3
import k3d

plot = k3d.plot(camera_auto_fit=False, grid_visible=False, 
                menu_visibility=True, axes_helper=False )

indices = mesh.simplices.astype(np.uint32)
points = np.column_stack(mesh.points.T).astype(np.float32)

mesh_viewer = k3d.mesh(points, indices, wireframe=False, attribute=data_n,
                   color_map=k3d.colormaps.paraview_color_maps.Cool_to_Warm, 
                   name="splines",
                   flat_shading=False, opacity=1.0  )


plot   += mesh_viewer
plot += background
plot.display()

## ## ## 

from ipywidgets import interact, interactive
import ipywidgets as widgets

choices = { "dlon": dlon1,
            "dlat": dlat1,
            "dlon + tension": dlon2,
            "dlat + tension": dlat2,
            "dlon delta": dlon1 - dlon2,
            "dlat delta": dlat1 - dlat2 }


@interact(choice=choices.keys())
def chooser(choice):
    mesh_viewer.attribute = choices[choice].astype(np.float32)
    range = choices[choice].max() * 0.2
    mesh_viewer.color_range = [-range, range]
    return 


```

The next notebook is [Ex9-Voronoi-Diagram](Ex9-Voronoi-Diagram.md)
