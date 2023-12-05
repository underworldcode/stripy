# Introduction to Stripy


[![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/underworldcode/stripy.svg)](https://hub.docker.com/r/underworldcode/stripy)
[![PyPI](https://img.shields.io/pypi/v/stripy.svg)](https://pypi.org/project/stripy/)

[![pip builds](https://github.com/underworldcode/stripy/workflows/pip%20builds/badge.svg)](https://github.com/underworldcode/stripy/actions?query=workflow%3A%22pip+builds%22)


A Python interface to TRIPACK and STRIPACK Fortran code for (constrained) triangulation in Cartesian coordinates and on a sphere. Stripy is an object-oriented package and includes routines from SRFPACK and SSRFPACK for interpolation (nearest neighbor, linear and hermite cubic) and to evaluate derivatives (Renka 1996a,b and 1997a,b).

---

<iframe src="images/k3d/icosohedron_example.html" title="Icosohedron" width="80%" height="500px"></iframe>

---

`stripy` is bundled with `litho1pt0` which is a python interface to the _crust 1.0_ dataset and the lithospheric part of the _litho 1.0_ dataset (Laske et al, 2013 and Pasyanos et al, 2014) which both requires and demonstrates the triangulation / searching and interpolation on the sphere that is provided by `stripy`.


```{figure} https://github.com/underworldcode/stripy/blob/master/stripy/Notebooks/Images/seafloor-age-topo.png?raw=true

---
width: 150mm
name: seafloor age mapped to equispaced points in the ocean
---
_Sample images created with `stripy` illustrating the meshing capability: ocean age data can be triangulated on the sphere with no need for points on land. Once `stripy` ingests your data points, you can sample another dataset to your grid (bathymetry on the right), smooth, find the derivatives of your data, or interpolate to another set of points._

```




## Links

  - [Stripy on Github](https://github.com/underworldcode/stripy)
  - [API documentation](https://underworldcode.github.io/stripy/2.0.5b2_api)



## Accessibility

<button type="button" onclick="legibleFontSwitcher()">Switch Font</button>&nbsp;&nbsp;<button type="button" onclick="fontScaler(1.1)">&#10133;</button><button type="button" onclick="fontScaler(0.0)">&#9679;</button><button type="button" onclick="fontScaler(0.909)">&#10134;</button>  


The online web page can also be typeset using the [Atkinson Hyperlegible](https://brailleinstitute.org/freefont) font everywhere, other than monospaced computer code, as an aid to legibility. This button is also located at the bottom of the left navigation bar on every page and will toggle between settings.

