## Interface to the Litho1.0 Earth model

Provides a simple, fast, interface to the litho1.0 model (Pasyanos et al, 2014) using the `stripack` spherical meshing / interpolation toolkit. It is a convenient replacement for the `access_litho` tool supplied with the data and considerably faster when computing values at lots of points.

This initial "release" provides layer depths for any (lon, lat) point on the Earth's surface, and depth-sampling of a particular quantity at a given (lon, lat) location.


## Litho 1.0

Pasyanos, M.E., T.G. Masters, G. Laske, and Z. Ma (2014). LITHO1.0: An updated crust and lithospheric model of the Earth, J. Geophys. Res., 119 (3), 2153-2173, DOI: 10.1002/2013JB010626..

Cover your eyes, then read the [litho 1.0 home page](http://igppweb.ucsd.edu/~gabi/litho1.0.html)

## Stripy

[Stripy](https://github.com/University-of-Melbourne-Geodynamics/stripy) is a wrapper to the FORTRAN libraries which carries out computational geometry tasks on the unit sphere, by Robert Renka. It is here used to provide fast location and interpolation through the Litho 1.0 dataset which is provided on a subdivided icosahedral spherical triangulation.
