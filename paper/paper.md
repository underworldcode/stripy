---
title: 'Stripy: A Python module for (constrained) triangulation in Cartesian coordinates and on a sphere.'
tags:
  - Python
  - triangulation
  - splines
  - interpolation
  - spherical mesh
  - Cartesian mesh
authors:
  - name: Ben Mather
    orcid: 0000-0003-3566-1557
    affiliation: 1
  - name: Louis Moresi
    orcid: 0000-0003-3685-174X
    affiliation: "2, 3" # (Multiple affiliations must be quoted)
affiliations:
 - name: School of Geoscience, The University of Sydney
   index: 1
 - name: School of Earth Science, The University of Melbourne
   index: 2
 - name: Research School of Earth Science, Australian National University
   index: 3
date: 7 March 2019
bibliography: paper.bib
---

![Figure](figure.png)

# Summary

The triangulation of scattered points is a commonly faced problem in all fields of science.
For most geographical applications, the spherical triangulation of unstructured points is required as most data are expected in longitude and latitude coordinates.
``Stripy`` is a Python interface to TRIPACK and STRIPACK Fortran code for (constrained) triangulation in Cartesian coordinates and on a sphere [@Renka1996a; @Renka1997a].
It is an object-oriented package and includes routines from SRFPACK and SSRFPACK for interpolation (nearest neighbour, linear, and hermite cubic) and to evaluate derivatives [@Renka1996b; @Renka1997b].
``Stripy`` offers a consistent API for interoperability between each coordinate system so users can transfer from Cartesian to spherical coordinates with minor changes to their code.

``Stripy`` includes the following functionality:

- Spherical and Cartesian triangulation of scattered points.
- Construction of Cartesian and Spherical meshes.
- Nearest-neighbour, linear, and hermite cubic interpolation.
- Evaluation of derivatives.
- Smoothing operations.
- Mesh refinement on line segments / triangle centroids.

These features are significant within Geographic Information Systems (GIS) where unstructured data can be triangulated, in projected or geographic coordinates, and interpolated to other datasets.
To that end, ``Stripy`` is bundled with ``litho1pt0``, a Python interface to the _crust 1.0_ dataset and the lithospheric part of the _litho 1.0_ dataset [@Laske2013; @Pasyanos2014], which both requires and demonstrates the triangulation, searching, and interpolation on the sphere that is provided by ``Stripy``.
The demonstrations operate in the Jupyter notebook environment for two matching sets of notebooks: one set for Cartesian triangulations, and one for spherical triangulations.

# Acknowledgements

Development of ``Stripy`` was financially supported by AuScope as part of the Simulation Analysis Modelling platform.

# References
