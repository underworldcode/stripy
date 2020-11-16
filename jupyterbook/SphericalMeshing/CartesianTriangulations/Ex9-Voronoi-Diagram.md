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

# Example 9 - Constructing the Voronoi diagram

The dual of a Delaunay triangulation is its Voronoi diagram. Stripy can efficiently calculate the voronoi points from a triangulation and, optionally, build the Voronoi regions for you.

## Contents

- [Voronoi points](#Voronoi-points)
- [Voronoi regions](#Voronoi-regions)

```{code-cell} ipython3
import numpy as np
import stripy
import matplotlib.pyplot as plt

%matplotlib inline
```

```{code-cell} ipython3
import stripy as stripy

xmin = 0.0
xmax = 10.0
ymin = 0.0
ymax = 10.0
extent = [xmin, xmax, ymin, ymax]

spacingX = 0.5
spacingY = 0.5

## Set up a coarse mesh for visual clarity

mesh = stripy.cartesian_meshes.elliptical_mesh(extent, spacingX, spacingY)
```

```{code-cell} ipython3
fig = plt.figure(figsize=(10, 8), facecolor="none")

ax = fig.add_subplot(111, title="Delaunay triangulation")
ax.axis('off')
ax.triplot(mesh.x, mesh.y, mesh.simplices, c='LightGrey', zorder=1)
ax.scatter(mesh.x, mesh.y, c='C0', marker='.', zorder=2)

plt.show()
```

## Voronoi points

A voronoi point (or circumcentre) exists for every triangle in the mesh. They are calculated by finding a constant radius (the circumradius) that is shared between each point in the triangle. Depending on the mesh, some voronoi points will reside outside the area contained within a triangle.

To calculate the Voronoi points from the mesh, use the `voronoi_points` method. This method may optionally return the circumradius, signed triangle area, or aspect ratio for each voronoi point.

```{code-cell} ipython3
vx, vy = mesh.voronoi_points()


fig = plt.figure(figsize=(10, 8), facecolor="none")

ax = fig.add_subplot(111, title="Delaunay triangulation + Voronoi points")
ax.axis('off')
ax.triplot(mesh.x, mesh.y, mesh.simplices, c='LightGrey', zorder=1)
ax.scatter(mesh.x, mesh.y, c='C0', marker='.', zorder=2)
ax.scatter(vx, vy, c='C1', marker='*', zorder=3)

plt.show()
```

## Voronoi regions

Often it is desirable to obtain the corresponding region enclosed by voronoi points. We can easily find the line segments connecting each voronoi point by retrieving the neighbour simplices.

```{code-cell} ipython3
neighbours = mesh.neighbour_simplices()

circumcentres = np.column_stack([vx, vy])
voronoi_edges = circumcentres[neighbours]
voronoi_edges[neighbours == -1] = np.nan # remove edges at infinity


lines = []
lines.extend(zip(circumcentres, voronoi_edges[:,0,:]))
lines.extend(zip(circumcentres, voronoi_edges[:,1,:]))
lines.extend(zip(circumcentres, voronoi_edges[:,2,:]))
```

```{code-cell} ipython3
# Plot it
from matplotlib.collections import LineCollection

linesC = LineCollection(lines, edgecolor='k')

fig = plt.figure(figsize=(10, 8), facecolor="none")
ax = fig.add_subplot(111, title='Voronoi diagram')
ax.axis('off')
ax.scatter(mesh.x, mesh.y, c='C0', marker='.', zorder=2)
ax.scatter(vx, vy, c='C1', marker='*', zorder=3)
ax.add_collection(linesC)

plt.show()
```

Looks good! But in order to find which voronoi points belong to a region, we need to iterate through all of the triangles. Since each vertex in the Delaunay represents a voronoi "site", we can accumulate all of the circumcentres into a list and sort them using `voronoi_points_and_regions`.

```{code-cell} ipython3
vx, vy, regions = mesh.voronoi_points_and_regions()


linesC = LineCollection(lines, edgecolor='k')

fig = plt.figure(figsize=(10, 8), facecolor="none")
ax = fig.add_subplot(111, title='Voronoi diagram')
ax.axis('off')
ax.scatter(mesh.x, mesh.y, c='C0', marker='.', zorder=2)
ax.scatter(vx, vy, c='C1', marker='*', zorder=3)
ax.add_collection(linesC)

# highlight specific region
r = 1
region = regions[r]

# iterate through region and connect up the points
for i in range(len(region)):
    i0 = region[i - 1]
    i1 = region[i]

    xx = [vx[i0], vx[i1]]
    yy = [vy[i0], vy[i1]]
    ax.plot(xx, yy, c='r', linewidth=5)
```
