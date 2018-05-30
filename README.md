karl
====

## Purpose

Karl, named after Karl Schwarzschild, is a library for numerical calculations regarding
black holes.

Current features:

* the Schwarzschild spacetime
* Kruskal-Szekeres coordinates 
* Schwarzschild coordinates
* Christoffel symbols
* test suite
* implementation in python
* computation of geodesics using Runge-Kutta integration

Partially complete:

* automatic switching between coordinate charts to avoid coordinate singularities in spherical coordinates and floating-point overflow in Kruskal-Szekeres coordinates
* graceful termination when a particle reaches a singularity

Features I plan to add:

* trajectories under the influence of a force (e.g., simulations of rocket ships)
* Reissner-Nordstrom spacetime
* spacetimes that lack spherical symmetry
* implementation in javascript

Documentation for the math is in the file doc.tex, which can be
compiled to pdf format by doing a "make doc." (Comments in the code do
not document the math or the definitions of the variables.)

## Installing

It's written in python3. Any changes required to make it work in python2 would probably be minimal,
but I haven't tried. It requires the scipy and numpy libraries, which can be installed on a debian
system by doing `apt-get install python-scipy3 python-numpy3`.

## Lack of test coverage

Test the code for transitioning a vector from Kruskal to Schwarzschild.

## Bugs

Currently the handle_transition method only knows about transitions
between Kruskal and Schwarzschild coordinates, need to take care of
the rot90 and era cases as well. 

