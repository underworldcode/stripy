Python interface to ACM algorithm 772 (STRIPACK) fortran code for triangulation
on a sphere.  Includes routines from algorithm 623 for interpolation
(nearest neighbor, linear and hermite cubic).

to install ([numpy](http://numpy.org) and fortran compiler (preferably 
[gfortran](https://gcc.gnu.org/wiki/GFortran)) required):

* ``python setup.py build``
   - If you change the fortran compiler, you may have to add the 
flags `config_fc --fcompiler=<compiler name>` when setup.py is run
(see docs for [numpy.distutils] (http://docs.scipy.org/doc/numpy-dev/f2py/distutils.html)).
* ``python setup.py install``

* to run test, execute ``python stripack/__init__.py``

* see source code <a href="http://htmlpreview.github.com/?https://github.com/jswhit/stripack/master/stripack.html" target="_blank">docstrings</a> for documentation, 
[``test/fib_test.py``](https://github.com/jswhit/stripack/blob/master/test/fib_test.py) for example usage.

References:

 R. J. Renka, "ALGORITHM 623:  Interpolation on the Surface of a
 Sphere", ACM Trans. Math. Software, Vol. 10, No. 4, December 1984,
 pp. 437-439. (http://dl.acm.org/citation.cfm?id=356107)

 R. J. Renka, "ALGORITHM 772: STRIPACK: Delaunay triangulation
 and Voronoi diagram on the surface of a sphere"
 ACM Trans. Math. Software, Volume 23 Issue 3, Sept. 1997
 pp 416-434. (http://dl.acm.org/citation.cfm?id=275329)
