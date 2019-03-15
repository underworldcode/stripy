# Stripy

[![Docker Cloud Automated build](https://img.shields.io/docker/cloud/automated/underworldcode/stripy.svg)](https://hub.docker.com/r/underworldcode/stripy)
[![PyPI](https://img.shields.io/pypi/v/stripy.svg)](https://pypi.org/project/stripy/)

A Python interface to TRIPACK and STRIPACK Fortran code for (constrained) triangulation in Cartesian coordinates and on a sphere. Stripy is an object-oriented package and includes routines from SRFPACK and SSRFPACK for interpolation (nearest neighbor, linear and hermite cubic) and to evaluate derivatives (Renka 1996a,b and 1997a,b).

`stripy` is bundled with `litho1pt0` which is a python interface to the _crust 1.0_ dataset and the lithospheric part of the _litho 1.0_ dataset (Laske et al, 2013 and Pasyanos et al, 2014) which both requires and demonstrates the triangulation / searching and interpolation on the sphere that is provided by `stripy`.


![Examples](https://github.com/underworldcode/stripy/blob/master/Notebooks/Images/Examples.png?raw=true)


_Sample images created with `stripy` illustrating the meshing capability, the ability to refine meshes to match criteria such as data density, and the ability to create distance-weighted averages to meshes and continuous interpolating functions_

#### Binder

Launch the demonstration at [mybinder.org](https://mybinder.org/v2/gh/underworldcode/stripy/binder?filepath=Notebooks%2F0-StartHere.ipynb)

[![badge](https://img.shields.io/badge/launch-stripy-E66581.svg?logo=data:image/png;base64,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)](https://mybinder.org/v2/gh/underworldcode/stripy/binder?filepath=Notebooks%2F0-StartHere.ipynb)


## Navigation / Notebooks


There are two matching sets of `stripy` notebooks - one set for [Cartesian Triangulations](#Cartesian) and one for [Spherical Triangulations](#Spherical). For most geographical applications, the spherical triangulations are the natural choice as they expect longitude and latitude coordinates (admittedly in radians).

Note: the Cartesian and Spherical notebooks can be obtained / installed from `stripy` itself as follows:

```bash
   python -c 'import stripy; stripy.documentation.install_documentation(path="Notebooks")'   
```

### Cartesian

  - [Ex1-Cartesian-Triangulations.ipynb](stripy-src/stripy/Notebooks/CartesianTriangulations/Ex1-Cartesian-Triangulations.ipynb)
  - [Ex2-CartesianGrids.ipynb](stripy-src/stripy/Notebooks/CartesianTriangulations/Ex2-CartesianGrids.ipynb)
  - [Ex3-Interpolation.ipynb](stripy-src/stripy/Notebooks/CartesianTriangulations/Ex3-Interpolation.ipynb)
  - [Ex4-Gradients.ipynb](stripy-src/stripy/Notebooks/CartesianTriangulations/Ex4-Gradients.ipynb)
  - [Ex5-Smoothing.ipynb](stripy-src/stripy/Notebooks/CartesianTriangulations/Ex5-Smoothing.ipynb)
  - [Ex6-Scattered-Data.ipynb](stripy-src/stripy/Notebooks/CartesianTriangulations/Ex6-Scattered-Data.ipynb)
  - [Ex7-Refinement-of-Triangulations.ipynb](stripy-src/stripy/Notebooks/CartesianTriangulations/Ex7-Refinement-of-Triangulations.ipynb)

### Spherical

  - [Ex1-Spherical-Triangulations.ipynb](stripy-src/stripy/Notebooks/SphericalTriangulations/Ex1-Spherical-Triangulations.ipynb)
  - [Ex2-SphericalGrids.ipynb](stripy-src/stripy/Notebooks/SphericalTriangulations/Ex2-SphericalGrids.ipynb)
  - [Ex3-Interpolation.ipynb](stripy-src/stripy/Notebooks/SphericalTriangulations/Ex3-Interpolation.ipynb)
  - [Ex4-Gradients.ipynb](stripy-src/stripy/Notebooks/SphericalTriangulations/Ex4-Gradients.ipynb)
  - [Ex5-Smoothing.ipynb](stripy-src/stripy/Notebooks/SphericalTriangulations/Ex5-Smoothing.ipynb)
  - [Ex6-Scattered-Data.ipynb](stripy-src/stripy/Notebooks/SphericalTriangulations/Ex6-Scattered-Data.ipynb)
  - [Ex7-Refinement-of-Triangulations.ipynb](stripy-src/stripy/Notebooks/SphericalTriangulations/Ex7-Refinement-of-Triangulations.ipynb)


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

To install ([numpy](http://numpy.org) and fortran compiler, preferably
[gfortran](https://gcc.gnu.org/wiki/GFortran), required):

```bash
python setup.py build
```
   - If you change the fortran compiler, you may have to add the
flags `config_fc --fcompiler=<compiler name>` when setup.py is run
(see docs for [numpy.distutils](http://docs.scipy.org/doc/numpy-dev/f2py/distutils.html)).
```bash
python setup.py install
```

Alternatively install using pip:

```bash
pip install [--user] stripy
```

## Usage

Two classes are included as part of the Stripy package:

- `sTriangulation` (Spherical coordinates)
- `Triangulation` (Cartesian coordinates)

These classes share similar methods and can be easily interchanged.
In addition, there are many helper functions provided for building meshes.

A series of tests are located in the *tests* subdirectory.


## Docker

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

## References

   1. Laske, G., G. Masters, and Z. Ma (2013), Update on CRUST1. 0—A 1-degree global model of Earth's crust, Geophys Research Abstracts, 15, EGU2013–2658.

   1. Pasyanos, M. E., T. G. Masters, G. Laske, and Z. Ma (2014), LITHO1.0: An updated crust and lithospheric model of the Earth, Journal of Geophysical Research-Solid Earth, 119(3), 2153–2173, doi:10.1002/2013JB010626.

   1. R. J. Renka, "ALGORITHM 751: TRIPACK: A Constrained Two- Dimensional Delaunay Triangulation Package" ACM Trans. Math. Software, Vol. 22, No. 1, 1996, pp. 1-8.

   1. R. J. Renka, "ALGORITHM 752: SRFPACK: Software for Scattered Data Fitting with a Constrained Surface under Tension", ACM Trans. Math. Software, Vol. 22, No. 1, 1996, pp. 9-17.

   1. R. J. Renka, "ALGORITHM 772: STRIPACK: Delaunay Triangulation and Voronoi Diagram on the Surface of a Sphere" ACM Trans. Math. Software, Vol. 23, No. 3, 1997, pp. 416-434.

   1. R. J. Renka, "ALGORITHM 773: SSRFPACK: Interpolation of Scattered Data on the Surface of a Sphere with a Surface under Tension", ACM Trans. Math. Software, Vol. 23, No. 3, 1997, pp. 437-439.
