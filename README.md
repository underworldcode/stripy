# Stripy

[![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/underworldcode/stripy.svg)](https://hub.docker.com/r/underworldcode/stripy)
[![PyPI](https://img.shields.io/pypi/v/stripy.svg)](https://pypi.org/project/stripy/)

[![pip builds](https://github.com/underworldcode/stripy/workflows/pip%20builds/badge.svg)](https://github.com/underworldcode/stripy/actions?query=workflow%3A%22pip+builds%22)

![Conda Deployment](https://github.com/underworldcode/stripy/workflows/Conda%20Deployment/badge.svg)


A Python interface to TRIPACK and STRIPACK Fortran code for (constrained) triangulation in Cartesian coordinates and on a sphere. Stripy is an object-oriented package and includes routines from SRFPACK and SSRFPACK for interpolation (nearest neighbor, linear and hermite cubic) and to evaluate derivatives (Renka 1996a,b and 1997a,b).

`stripy` is bundled with `litho1pt0` which is a python interface to the _crust 1.0_ dataset and the lithospheric part of the _litho 1.0_ dataset (Laske et al, 2013 and Pasyanos et al, 2014) which both requires and demonstrates the triangulation / searching and interpolation on the sphere that is provided by `stripy`.


![Examples](https://github.com/underworldcode/stripy/blob/master/stripy/Notebooks/Images/seafloor-age-topo.png?raw=true)

_Sample images created with `stripy` illustrating the meshing capability: ocean age data can be triangulated on the sphere with no need for points on land. Once `stripy` ingests your data points, you can sample another dataset to your grid (bathymetry on the right), smooth, find the derivatives of your data, or interpolate to another set of points._

## Documentation

There are two matching sets of `stripy` notebooks - one set for [Cartesian Triangulations](#Cartesian) and one for [Spherical Triangulations](#Spherical). For most geographical applications, the spherical triangulations are the natural choice as they expect longitude and latitude coordinates (admittedly in radians).  There are some worked examples
which use the companion package litho1pt0

### Stable code 

  - Documentation / Notebooks [https://underworldcode.github.io/stripy/2.0.5b2](https://underworldcode.github.io/stripy/2.0.5b2)
  - API documentation [https://underworldcode.github.io/stripy/2.0.5b2_api](https://underworldcode.github.io/stripy/2.0.5b2_api)

### Bleeding edge code 


  - Documentation / Notebooks [https://underworldcode.github.io/stripy/2.1.0b1](https://underworldcode.github.io/stripy/2.1.0b1)
  - API documentation [https://underworldcode.github.io/stripy/2.1.0b1_api](https://underworldcode.github.io/stripy/2.1.0b1_api)

For previous versions, see the [changelog](Changelog.md)

### Installation & Running in the cloud

#### Binder

Launch the demonstration <!-- at [links.underworldcode.org/stripy-live (mybinder.org)](http://links.underworldcode.org/stripy-live) -->

[![Binder](https://mybinder.org/badge_logo.svg)](http://links.underworldcode.org/stripy-live)

(This is the development branch)

[![Binder-dev](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/underworldcode/stripy/dev?filepath=Notebooks%2F0-StartHere.ipynb)


#### Citation

[![DOI](http://joss.theoj.org/papers/10.21105/joss.01410/status.svg)](https://doi.org/10.21105/joss.01410)

Moresi, L. and Mather, B.R., (2019). Stripy: A Python module for (constrained) triangulation in Cartesian coordinates and on a sphere.. Journal of Open Source Software, 4(38), 1410, [https://doi.org/10.21105/joss.01410](https://doi.org/10.21105/joss.01410)

## Navigation / Notebooks


There are two matching sets of `stripy` notebooks - one set for [Cartesian Triangulations](#Cartesian) and one for [Spherical Triangulations](#Spherical). For most geographical applications, the spherical triangulations are the natural choice as they expect longitude and latitude coordinates (admittedly in radians).

Note: the Cartesian and Spherical notebooks can be obtained / installed from `stripy` itself as follows:

```bash
   python -c 'import stripy; stripy.documentation.install_documentation(path="Notebooks")'
```

### Cartesian

  - [Ex1-Cartesian-Triangulations.ipynb](stripy/Notebooks/CartesianTriangulations/Ex1-Cartesian-Triangulations.ipynb)
  - [Ex2-CartesianGrids.ipynb](stripy/Notebooks/CartesianTriangulations/Ex2-CartesianGrids.ipynb)
  - [Ex3-Interpolation.ipynb](stripy/Notebooks/CartesianTriangulations/Ex3-Interpolation.ipynb)
  - [Ex4-Gradients.ipynb](stripy/Notebooks/CartesianTriangulations/Ex4-Gradients.ipynb)
  - [Ex5-Smoothing.ipynb](stripy/Notebooks/CartesianTriangulations/Ex5-Smoothing.ipynb)
  - [Ex6-Scattered-Data.ipynb](stripy/Notebooks/CartesianTriangulations/Ex6-Scattered-Data.ipynb)
  - [Ex7-Refinement-of-Triangulations.ipynb](stripy/Notebooks/CartesianTriangulations/Ex7-Refinement-of-Triangulations.ipynb)
  - [Ex8-Spline-Tension.ipynb](stripy/Notebooks/CartesianTriangulations/Ex8-Spline-Tension.ipynb)
  - [Ex9-Voronoi-Diagram.ipynb](stripy/Notebooks/CartesianTriangulations/Ex9-Voronoi-Diagram.ipynb)

### Spherical

  - [Ex1-Spherical-Triangulations.ipynb](stripy/Notebooks/SphericalTriangulations/Ex1-Spherical-Triangulations.ipynb)
  - [Ex2-SphericalGrids.ipynb](stripy/Notebooks/SphericalTriangulations/Ex2-SphericalGrids.ipynb)
  - [Ex3-Interpolation.ipynb](stripy/Notebooks/SphericalTriangulations/Ex3-Interpolation.ipynb)
  - [Ex4-Gradients.ipynb](stripy/Notebooks/SphericalTriangulations/Ex4-Gradients.ipynb)
  - [Ex5-Smoothing.ipynb](stripy/Notebooks/SphericalTriangulations/Ex5-Smoothing.ipynb)
  - [Ex6-Scattered-Data.ipynb](stripy/Notebooks/SphericalTriangulations/Ex6-Scattered-Data.ipynb)
  - [Ex7-Refinement-of-Triangulations.ipynb](stripy/Notebooks/SphericalTriangulations/Ex7-Refinement-of-Triangulations.ipynb)
  - [Ex8-Spline-Tension.ipynb](stripy/Notebooks/SphericalTriangulations/Ex8-Spline-Tension.ipynb)
  - [Ex9-Voronoi-Diagram.ipynb](stripy/Notebooks/SphericalTriangulations/Ex9-Voronoi-Diagram.ipynb)


### Examples

Note, these examples are the notebooks from `litho1pt0` which are installed from the
package itself:

```bash
   python -c 'import litho1pt0; litho1pt0.documentation.install_documentation(path="Notebooks")'
```

The first three notebooks are an introduction to `litho1pt0` that does not explicitly mention `stripy` but
the next two worked examples show how to search, interpolate and plot with the help of `stripy` routines.

  - [Ex1-Litho1Layers.ipynb](https://github.com/underworldcode/litho1pt0/blob/master/litho1pt0/Notebooks/litho1pt0/Ex1-Litho1Layers.ipynb)
  - [Ex2-Litho1Properties.ipynb](https://github.com/underworldcode/litho1pt0/blob/master/litho1pt0/Notebooks/litho1pt0/Ex2-Litho1Properties.ipynb)
  - [Ex3-CrustalRegionalisation.ipynb](https://github.com/underworldcode/litho1pt0/blob/master/litho1pt0/Notebooks/litho1pt0/Ex3-CrustalRegionalisation.ipynb)
  - [WorkEx1-CratonAverageProperties.ipynb](https://github.com/underworldcode/litho1pt0/blob/master/litho1pt0/Notebooks/litho1pt0/WorkEx1-CratonAverageProperties.ipynb)
  - [WorkEx2-OceanDepthAge.ipynb](https://github.com/underworldcode/litho1pt0/blob/master/litho1pt0/Notebooks/litho1pt0/WorkEx2-OceanDepthAge.ipynb)


## Installation

### Dependencies

You will need **Python 3.6+**.
Also, the following packages are required:

 - [`gfortran`](https://www.fatiando.org/verde/latest/install.html) (or any Fortran compiler)
 - [`numpy`](http://numpy.org)
 - [`scipy`](https://scipy.org)

**Recommended Packages** for running the notebooks:

 - [`litho1pt0`](https://pypi.org/project/litho1pt0/)
 - [`matplotlib`](https://matplotlib.org/)
 - [`imageio`](https://imageio.github.io/)
 - [`cartopy`](https://scitools.org.uk/cartopy/docs/latest/)
 - [`k3d`](https://github.com/K3D-tools/K3D-jupyter)
 - [`xarray`](http://xarray.pydata.org/en/stable/)
 - [`netcdf4`](https://unidata.github.io/netcdf4-python/)

All of which should be available from pip or anaconda (conda-forge) for most platforms.

### Installing using pip

You can install `stripy` using the
[`pip package manager`](https://pypi.org/project/pip/) with either version of Python:

```bash
python3 -m pip install stripy
```

All the dependencies will be automatically installed by `pip`, except for `gfortran`
(or any Fortran compiler). It must be installed in your system before installing
`stripy` with `pip`.

If you change the Fortran compiler, you may have to add the
flags `config_fc --fcompiler=<compiler name>` when `setup.py` is run
(see docs for [numpy.distutils](http://docs.scipy.org/doc/numpy-dev/f2py/distutils.html)).

### Installing with conda

If you use the anaconda packaging system, then you should be able to 

```bash
conda install -c geo-down-under stripy
```

## Usage

Two classes are included as part of the Stripy package:

- `sTriangulation` (Spherical coordinates)
- `Triangulation` (Cartesian coordinates)

These classes share similar methods and can be easily interchanged.
In addition, there are many helper functions provided for building meshes.

A series of tests are located in the *tests* subdirectory.
In order to perform these tests clone the repository and run [`pytest`](https://pypi.org/project/pytest/):

```bash
git checkout https://github.com/underworldcode/stripy.git
cd stripy
pytest -v
```


## References

   1. Laske, G., G. Masters, and Z. Ma (2013), Update on CRUST1. 0—A 1-degree global model of Earth's crust, Geophys Research Abstracts, 15, EGU2013–2658.

   1. Pasyanos, M. E., T. G. Masters, G. Laske, and Z. Ma (2014), LITHO1.0: An updated crust and lithospheric model of the Earth, Journal of Geophysical Research-Solid Earth, 119(3), 2153–2173, doi:10.1002/2013JB010626.

   1. R. J. Renka, "ALGORITHM 751: TRIPACK: A Constrained Two- Dimensional Delaunay Triangulation Package" ACM Trans. Math. Software, Vol. 22, No. 1, 1996, pp. 1-8.

   1. R. J. Renka, "ALGORITHM 752: SRFPACK: Software for Scattered Data Fitting with a Constrained Surface under Tension", ACM Trans. Math. Software, Vol. 22, No. 1, 1996, pp. 9-17.

   1. R. J. Renka, "ALGORITHM 772: STRIPACK: Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere" ACM Trans. Math. Software, Vol. 23, No. 3, 1997, pp. 416-434.

   1. R. J. Renka, "ALGORITHM 773: SSRFPACK: Interpolation of Scattered Data on the Surface of a Sphere with a Surface under Tension", ACM Trans. Math. Software, Vol. 23, No. 3, 1997, pp. 437-439.
