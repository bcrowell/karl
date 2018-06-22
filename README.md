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
* implementation in both python and javascript
* computation of geodesics using Runge-Kutta integration

Features I plan to add:

* automatic switching between coordinate charts for better numerical behavior
* automatic handling of the case where a geodesic is terminating
* trajectories under the influence of a force (e.g., simulations of rocket ships)
* Reissner-Nordstrom spacetime

Documentation for the math is in the file doc.tex, which can be
compiled to pdf format by doing a "make doc." (Comments in the code do
not document the math or the definitions of the variables.)

## Installing

It's written in python3 and requires the scipy and numpy libraries, and the filepp utility,
which can be installed on a debian
system by doing `apt-get install python-scipy3 python-numpy3 filepp`.
To translate python to javascript:

    apt-get install rhino
    pip install jsbeautifier
