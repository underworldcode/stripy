# Stripy

A Python interface to TRIPACK and STRIPACK Fortran code for (constrained) triangulation in Cartesian coordinates and on a sphere. Stripy is an object-oriented package and includes routines from SRFPACK and SSRFPACK for interpolation (nearest neighbor, linear and hermite cubic) and to evaluate derivatives.

This package was inspired by Jeffrey Whittaker's [implementation](https://github.com/jswhit/stripack) for STRIPACK.

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

- sTriangulation (Spherical coordinates)
- Triangulation (Cartesian coordinates)

These classes share similar methods and can be easily interchanged.
In addition, there are many helper functions provided for building

A series of tests are located in the *tests* subdirectory.

## Dockerfile

To build the dockerfile run the following from the home stripy directory

```bash
source ./Docker/scripts/build-base-dockerfile.sh`
```

Then either run this from the kitematic gui or using

```bash
docker run -P -i lmoresi/stripy:0.30
```

This will serve the jupyter notebook browser on a free port.
Kitematic should provide a preview that can be clicked to open the demonstration notebooks in a browser.


## References

 R. J. Renka, "ALGORITHM 751: TRIPACK: A Constrained Two-
 Dimensional Delaunay Triangulation Package" ACM Trans. Math.
 Software, Vol. 22, No. 1, 1996, pp. 1-8.

 R. J. Renka, "ALGORITHM 752: SRFPACK: Software for Scattered
 Data Fitting with a Constrained Surface under Tension", ACM
 Trans. Math. Software, Vol. 22, No. 1, 1996, pp. 9-17.

 R. J. Renka, "ALGORITHM 772: STRIPACK: Delaunay Triangulation
 and Voronoi Diagram on the Surface of a Sphere"
 ACM Trans. Math. Software, Vol. 23, No. 3, 1997, pp. 416-434.

 R. J. Renka, "ALGORITHM 773: SSRFPACK: Interpolation of Scattered
 Data on the Surface of a Sphere with a Surface under Tension",
 ACM Trans. Math. Software, Vol. 23, No. 3, 1997, pp. 437-439.
