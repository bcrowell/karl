#!/usr/bin/python3

#include "math.h"
#include "init.h"
#include "test.h"

import schwarzschild

# Exercise as many lines of code as possible, using random coordinates.
t = 3.3
r = 2.2
sin_theta = 1/sqrt(2.0)
cos_theta = 1/sqrt(2.0)
# A set of i,j,k that lie on the unit sphere:
i=0.6
j=0.8
k=0.0

schwarzschild.metric_sch4(r,sin_theta)
schwarzschild.sigma(r)
schwarzschild.christoffel_sch4(t,r,sin_theta,cos_theta)
schwarzschild.christoffel([t,r,i,j,k])

# As far as testing whether the results are actually correct, see
# tests in test_runge_kutta.

done(verbosity,"schwarzschild")
