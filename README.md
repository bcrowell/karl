karl
====

## Purpose

Karl, named after Karl Schwarzschild, is a library for numerical calculations of trajectories
of test particles in black-hole spacetimes. The technique used is direct numerical integration
of the equations of motion. Compared to techniques involving conserved quantities, this has
some advantages:

* We can simulate a nongravitational force, e.g., in the motion of a rocket ship.
* We can handle spacetimes having fewer symmetries than the Schwarzschild spacetime.
* We can do calculations involving proper time or other affine parameters, e.g., finding how much time it takes for a particle to hit the singularity.

Current features:

* the Schwarzschild spacetime
* computation of trajectories using Runge-Kutta integration
* Schwarzschild coordinates
* Kruskal-Szekeres coordinates 
* Christoffel symbols, coordinate transformations, and Jacobian matrices
* test suite
* implementation in both python and javascript, with time-critical code written in C for speed in the python version
* user-defined triggers to stop the computation when a specified coordinate or velocity reaches a set value

Features I plan to add:

* automatic switching between coordinate charts for better numerical behavior
* automatic handling of the case where a geodesic is terminating
* other black-hole spacetimes such as Reissner-Nordstrom

Documentation for the math is in the file doc.tex, which can be
compiled to pdf format by doing a "make doc." (Comments in the code do
not document the math or the definitions of the variables.)

## Installing

It's written in python3 and requires the scipy and numpy libraries, and the filepp utility,
which can be installed on a debian
system by doing `apt-get install python-scipy3 python-numpy3 filepp`.
If you want to be able to modify the python code and then translate python to javascript:

    apt-get install rhino ruby
    pip install jsbeautifier

## In the browser

There is a simple text-based demo in the file browser/sample.html.

The startup scripts in browser/util create the following global variables:
karl, print, verbosity, and some all-caps constants in constants.js.

## Notes on performance

Surprisingly (for me), the javascript implementation is an order of magnitude faster than
the python version, even after the time-critical python code was rewritten in C for speed.

Typical users of the js version will probably be using it in a browser, but if using it
on the command line through rhino, an extra speed boost can be achieved by running
rhino with the option -opt 9. This is not done by the tests by default, because it
prevents debugging information from being printed.

## Bugs

In test_kruskal, the test test_motion_kruskal_vs_schwarzschild fails due to poor
precision if we replace n=100 with n=1000. Is this 4th-order behavior with a larger
constant of proportionality than I'd imagined, or is it a bug or a sign of numerical
instability in Kruskal coordinates?

## To do

Wrap the simple Runge-Kutta routine in a fancier one that is adaptive, so that we can,
e.g., accurately determine the affine parameter at which we hit the singularity.
Also, the adaptive routine should automatically switch coordinate systems when
appropriate. Kruskal coordinates should be used only near the horizon.
(Disadvantages of Kruskal coordinates: computationally expensive, computationally
expensive in arcsinh-Kruskal to determine when we've hit the singularity.)

Allow sanity tests for whether the geodesic maintains its timelike, null, or
spacelike character, and in the timelike case for whether it stays in the future
or past light cone. A common numerical behavior when hitting the singularity is to get a
bounce.

## Ideas for samples and educational apps

Animation of a cloud of test particles undergoing tidal distortions.

Animation on a Penrose diagram.
Show a spacelike geodesic crossing from region I to region III, or timelike geodesics
from I and III meeting in II.

Make a 2-d spacewars-style game.

## Adding a new coordinate system

Define a number for it, and add that to spacetimes.h and spacetimes_c.h.
Calculate the Christoffel symbols as in, e.g., kruskal5.mac, and run the output
through clean_up_christoffel.rb. Implement the christoffel symbols in python.
Add the new python file to gen_depends.py.
Add to runge_kutta.chart_info.
Add to transform.chart_info() and,
if necessary,  to other routines in transform module.
