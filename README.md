# Stripy

[![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/underworldcode/stripy.svg)](https://hub.docker.com/r/underworldcode/stripy)
[![PyPI](https://img.shields.io/pypi/v/stripy.svg)](https://pypi.org/project/stripy/)
[![Build Status (master)](https://travis-ci.org/underworldcode/stripy.svg?branch=master)](https://travis-ci.org/underworldcode/stripy)

A Python interface to TRIPACK and STRIPACK Fortran code for (constrained) triangulation in Cartesian coordinates and on a sphere. Stripy is an object-oriented package and includes routines from SRFPACK and SSRFPACK for interpolation (nearest neighbor, linear and hermite cubic) and to evaluate derivatives (Renka 1996a,b and 1997a,b).

`stripy` is bundled with `litho1pt0` which is a python interface to the _crust 1.0_ dataset and the lithospheric part of the _litho 1.0_ dataset (Laske et al, 2013 and Pasyanos et al, 2014) which both requires and demonstrates the triangulation / searching and interpolation on the sphere that is provided by `stripy`.


![Examples](https://github.com/underworldcode/stripy/blob/master/stripy/Notebooks/Images/seafloor-age-topo.png?raw=true)

_Sample images created with `stripy` illustrating the meshing capability: ocean age data can be triangulated on the sphere with no need for points on land. Once `stripy` ingests your data points, you can sample another dataset to your grid (bathymetry on the right), smooth, find the derivatives of your data, or interpolate to another set of points._


#### Binder

Launch the demonstration at [links.underworldcode.org/stripy-live (mybinder.org)](http://links.underworldcode.org/stripy-live)

[![Binder](https://mybinder.org/badge_logo.svg)](http://links.underworldcode.org/stripy-live)

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

### Spherical

  - [Ex1-Spherical-Triangulations.ipynb](stripy/Notebooks/SphericalTriangulations/Ex1-Spherical-Triangulations.ipynb)
  - [Ex2-SphericalGrids.ipynb](stripy/Notebooks/SphericalTriangulations/Ex2-SphericalGrids.ipynb)
  - [Ex3-Interpolation.ipynb](stripy/Notebooks/SphericalTriangulations/Ex3-Interpolation.ipynb)
  - [Ex4-Gradients.ipynb](stripy/Notebooks/SphericalTriangulations/Ex4-Gradients.ipynb)
  - [Ex5-Smoothing.ipynb](stripy/Notebooks/SphericalTriangulations/Ex5-Smoothing.ipynb)
  - [Ex6-Scattered-Data.ipynb](stripy/Notebooks/SphericalTriangulations/Ex6-Scattered-Data.ipynb)
  - [Ex7-Refinement-of-Triangulations.ipynb](stripy/Notebooks/SphericalTriangulations/Ex7-Refinement-of-Triangulations.ipynb)


### Examples

Note, these examples are the notebooks from `litho1pt0` which are installed from the
package itself:

```bash
   python -c 'import litho1pt0; litho1pt0.documentation.install_documentation(path="Notebooks")'
```

The first three notebooks are an introduction to `litho1pt0` that does not explicitly mention `stripy` but
the next two worked examples show how to search, interpolate and plot with the help of `stripy` routines.

  - [Ex1-Litho1Layers.ipynb](litho1pt0-src/litho1pt0/Notebooks/litho1pt0/Ex1-Litho1Layers.ipynb)
  - [Ex2-Litho1Properties.ipynb](litho1pt0-src/litho1pt0/Notebooks/litho1pt0/Ex2-Litho1Properties.ipynb)
  - [Ex3-CrustalRegionalisation.ipynb](litho1pt0-src/litho1pt0/Notebooks/litho1pt0/Ex3-CrustalRegionalisation.ipynb)
  - [WorkEx1-CratonAverageProperties.ipynb](litho1pt0-src/litho1pt0/Notebooks/litho1pt0/WorkEx1-CratonAverageProperties.ipynb)
  - [WorkEx2-OceanDepthAge.ipynb](litho1pt0-src/litho1pt0/Notebooks/litho1pt0/WorkEx2-OceanDepthAge.ipynb)


## Installation

### Dependencies

You will need **Python 2.7 or 3.5+**.
Also, the following packages are required:

 - [`gfortran`](https://www.fatiando.org/verde/latest/install.html) (or any Fortran compiler)
 - [`numpy`](http://numpy.org)
 - [`scipy`](https://scipy.org)

**Recommended Packages** for running the notebooks:

 - [`litho1pt0`](https://pypi.org/project/litho1pt0/)
 - [`gdal`](https://www.gdal.org/)
 - [`matplotlib`](https://matplotlib.org/)
 - [`imageio`](https://imageio.github.io/)
 - [`cartopy`](https://scitools.org.uk/cartopy/docs/latest/)
 - [`pyproj`](https://github.com/pyproj4/pyproj)
 - [`lavavu`](https://github.com/OKaluza/LavaVu/)

### Installing using pip

You can install `stripy` using the
[`pip package manager`](https://pypi.org/project/pip/) with either version of Python:

```bash
python2 -m pip install stripy
python3 -m pip install stripy
```

All the dependencies will be automatically installed by `pip`, except for `gfortran`
(or any Fortran compiler). It must be installed in your system before installing
`stripy` with `pip`.

If you change the Fortran compiler, you may have to add the
flags `config_fc --fcompiler=<compiler name>` when `setup.py` is run
(see docs for [numpy.distutils](http://docs.scipy.org/doc/numpy-dev/f2py/distutils.html)).

### Installing using Docker

A more straightforward installation which does not depend on specific compilers relies on the [docker](http://www.docker.com) virtualisation system.

To install the docker image and test it is working:

```bash
   docker pull underworldcode/stripy:latest
   docker run --rm underworldcode/stripy:latest help
```

To install the helper scripts for bash:

```bash
   docker run --rm underworldcode/stripy:latest bash_utils > bash_utils.sh
   source bash_utils.sh
```

( you may find it helpful to move/rename this file and source it from
  your bash profile at login time )

The bash_utils.sh script installs the following functions which are
available through the bash command line:

```bash
  stripy-docker-help
  stripy-docker-sh
  stripy-docker-nb
  stripy-docker-browse
  stripy-docker-serve
  stripy-docker-terminal
```

For more information on these functions, run

```bash
  source bash_utils.sh
  stripy-docker-help
```

To use the docker version as you would, say, using ipython to type on the command line:

```bash
   source bash_utils.sh  # (only needs to be done once)
   stripy-docker-terminal
   ls
   ipython
```

To use the docker version to run a script

```bash
   source bash_utils.sh  # (only needs to be done once)
   stripy-docker-sh my_python_script.py
```

To build the dockerfile locally, we provide a script. First ensure you have checked out the source code from github and then run the script in the Docker directory. If you modify the dockerfile and want to push the image to make it publicly available, it will need to be retagged to upload somewhere other than the underworldcode repository.

```bash
git checkout https://github.com/underworldcode/stripy.git
cd stripy
source Docker/build-dockerfile.sh
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

### API Documentation

The API for all functions and classes in `stripy` can be accessed from [https://underworldcode.github.io/stripy/](https://underworldcode.github.io/stripy/).


## References

   1. Laske, G., G. Masters, and Z. Ma (2013), Update on CRUST1. 0—A 1-degree global model of Earth's crust, Geophys Research Abstracts, 15, EGU2013–2658.

   1. Pasyanos, M. E., T. G. Masters, G. Laske, and Z. Ma (2014), LITHO1.0: An updated crust and lithospheric model of the Earth, Journal of Geophysical Research-Solid Earth, 119(3), 2153–2173, doi:10.1002/2013JB010626.

   1. R. J. Renka, "ALGORITHM 751: TRIPACK: A Constrained Two- Dimensional Delaunay Triangulation Package" ACM Trans. Math. Software, Vol. 22, No. 1, 1996, pp. 1-8.

   1. R. J. Renka, "ALGORITHM 752: SRFPACK: Software for Scattered Data Fitting with a Constrained Surface under Tension", ACM Trans. Math. Software, Vol. 22, No. 1, 1996, pp. 9-17.

   1. R. J. Renka, "ALGORITHM 772: STRIPACK: Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere" ACM Trans. Math. Software, Vol. 23, No. 3, 1997, pp. 416-434.

   1. R. J. Renka, "ALGORITHM 773: SSRFPACK: Interpolation of Scattered Data on the Surface of a Sphere with a Surface under Tension", ACM Trans. Math. Software, Vol. 23, No. 3, 1997, pp. 437-439.
