karl
=====

Karl, named after Karl Schwarzschild, is a library for numerical calculations regarding
black holes.

Current features:

* the Schwarzschild spacetime
* Kruskal-Szekeres coordinates 
* Schwarzschild coordinates
* automatic switching between coordinate charts to avoid coordinate singularities in spherical coordinates and floating-point overflow in Kruskal-Szekeres coordinates
* Christoffel symbols
* test suite
* implementation in python

Features I plan to add:

* Reissner-Nordstrom spacetime
* spacetimes that lack spherical symmetry
* implementation in javascript

Documentation for the math is in the file doc.tex, which can be
compiled to pdf format by doing a "make doc." (Comments in the code do
not document the math or the definitions of the variables.)


